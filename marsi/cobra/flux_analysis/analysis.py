# Copyright 2016 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import logging

from IProgress import ProgressBar, Bar, Percentage
from bokeh.layouts import column
from bokeh.models import FactorRange, Range1d, LinearAxis
from bokeh.plotting import figure, show
from bokeh.charts import Line
from cobra.core import Reaction
from numpy import array

from pandas import DataFrame

from cobra.core import Metabolite
from cameo.core.solver_based_model import SolverBasedModel as Model
from cameo.core.result import Result
from cameo.util import TimeMachine
from cameo.flux_analysis.analysis import flux_variability_analysis, FluxVariabilityResult
from cameo.flux_analysis.simulation import pfba, fba
from cameo.exceptions import SolveError, Infeasible

from marsi.cobra.flux_analysis.manipulation import knockout_metabolite, compete_metabolite, inhibit_metabolite
from marsi.utils import frange


__all__ = ['metabolite_knockout_fitness', 'metabolite_knockout_phenotype', 'sensitivity_analysis']


BASE_ELEMENTS = ["C", "N"]


logger = logging.getLogger(__name__)


class MetaboliteKnockoutFitness(Result):
    """
    The result holder for the fitness landscape analysis.

    See Also
    --------
    marsi.cobra.flux_analysis.analysis.metabolite_knockout_fitness

    Attributes
    ----------
    data_frame : pandas.DataFrame
        A DataFrame with the result.


    """
    def __init__(self, fitness_data_frame, *args, **kwargs):
        super(MetaboliteKnockoutFitness, self).__init__(*args, **kwargs)
        assert isinstance(fitness_data_frame, DataFrame)
        self._data_frame = fitness_data_frame

    @property
    def data_frame(self):
        return DataFrame(self._data_frame)

    def plot(self, grid=None, width=None, height=None, title=None, *args, **kwargs):
        """
        """
        data = self._data_frame.sort_values('fitness')
        data['x'] = data.index
        show(Line(data, 'x', 'fitness', title=title, plot_width=width, plot_height=height))

    def _repr_html_(self):
        return self.plot(height=500, width=12 * len(self._data_frame), title="Metabolite knockout fitness landscape")


def metabolite_knockout_fitness(model, simulation_method=pfba, compartments=None, elements=BASE_ELEMENTS,
                                objective=None, ndecimals=6, progress=False, ncarbons=2, **simulation_kwargs):
    """
    Calculate the landscape of fitness for each metabolite knockout in the model.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        A constraint-based model.
    simulation_method : cameo.flux_analysis.simulation.fba
        A method to simulate the knockouts (e.g. cameo.flux_analysis.simulation.moma)
    compartments : list
        The compartments to consider (e.g. ["c", "g", "r"] for cytosol, golgi aparatus and endoplasmatic reticulum)
    elements : list
        Atomic elements to add to the result.
    objective : str, cameo.core.Reaction, other
        A valid objective for the model.
    ndecimals : int
        Number of decimals to use as precision.
    progress : bool
        Report progress.
    ncarbons : int
        Minimum number of carbons to consider.
    simulation_kwargs : dict
        Arguments for `simulation_method

    Returns
    -------
    MetaboliteKnockoutFitness
        The fitness landscape.
    """
    assert isinstance(model, Model)
    fitness = DataFrame(columns=["fitness"] + list(elements))

    if compartments is None:
        compartments = list(model.compartments.keys())

    if progress:
        iterator = ProgressBar(maxval=len(model.metabolites), widgets=[Bar(), Percentage()])
    else:
        iterator = iter
    for met in iterator(model.metabolites):
        if met.compartment in compartments and met.elements.get("C", 0) > ncarbons:
            with TimeMachine() as tm:
                knockout_metabolite(model, met, allow_accumulation=True, ignore_transport=True, time_machine=tm)
                try:
                    solution = simulation_method(model, objective=objective, **simulation_kwargs)
                    fitness.loc[met.id] = [round(solution[objective], ndecimals)] + \
                                          [met.elements.get(el, 0) for el in elements]
                except SolveError:
                    fitness.loc[met.id] = [.0] + [met.elements.get(el, 0) for el in elements]

    return MetaboliteKnockoutFitness(fitness)


class MetaboliteKnockoutPhenotypeResult(MetaboliteKnockoutFitness):
    def __init__(self, phenotype_data_frame, *args, **kwargs):
        super(MetaboliteKnockoutPhenotypeResult, self).__init__(phenotype_data_frame, *args, **kwargs)

    def __getitem__(self, item):
        phenotype = {}
        if isinstance(item, (int, slice)):
            fva = self._data_frame['fva'].iloc[item]
        elif isinstance(item, str):
            fva = self._data_frame['fva'].loc[item]
        elif isinstance(item, Metabolite):
            fva = self._data_frame['fva'].loc[item.id]
        else:
            raise ValueError("%s is not a valid search retrieval")

        for reaction_id, row in fva.data_frame.iterrows():
            if row['upper_bound'] == 0 and row['lower_bound'] == 0:
                continue
            phenotype[reaction_id] = (row['upper_bound'], row['lower_bound'])

        return phenotype

    def _phenotype_plot(self, index, factors):
        plot = figure(title=index, y_range=FactorRange(factors))
        phenotype = self[index]

        y0 = [phenotype.get(f, [0, 0])[0] for f in factors]
        y1 = [phenotype.get(f, [0, 0])[1] for f in factors]

        x0 = [(i - .5) for i in range(len(factors))]
        x1 = [(i + .5) for i in range(len(factors))]

        plot.quad(x0, x1, y0, y1)

        return plot

    def plot(self, indexes=None, grid=None, width=None, height=None, title=None, conditions=None, *args, **kwargs):
        data = self.data_frame
        if conditions:
            data = data.query(conditions)
        if indexes:
            data = data.loc[indexes]

        factors = list(set(sum(data.phenotype.apply(lambda p: list(p.keys())), [])))
        plots = []

        for index in data.index:
            plots.append(self._phenotype_plot(index, factors))
        show(column(plots))


def metabolite_knockout_phenotype(model, compartments=None, objective=None, ndecimals=6, elements=BASE_ELEMENTS,
                                  progress=False, ncarbons=2):
    assert isinstance(model, Model)
    phenotype = DataFrame(columns=['fitness', 'fva'] + elements)
    exchanges = model.exchanges

    if progress:
        iterator = ProgressBar(maxval=len(model.metabolites), widgets=[Bar(), Percentage()])
    else:
        iterator = iter

    for met in iterator(model.metabolites):
        if met.compartment in compartments and met.elements.get("C", 0) > ncarbons:
            with TimeMachine() as tm:
                knockout_metabolite(model, met, allow_accumulation=True, ignore_transport=True, time_machine=tm)
                fitness = fba(model, objective=objective)
                fva = flux_variability_analysis(model, reactions=exchanges, fraction_of_optimum=1)
                fva = FluxVariabilityResult(fva.data_frame.apply(round, args=(ndecimals,)))
                phenotype.loc[met.id] = [fitness, fva] + [met.elements.get(el, 0) for el in elements]

    return MetaboliteKnockoutPhenotypeResult(phenotype)


class SensitivityAnalysisResult(Result):
    def __init__(self, species_id, exchange_fluxes, steps, is_essential, biomass_fluxes=None,
                 biomass=None, variables=None, variables_fluxes=None, *args, **kwargs):
        super(SensitivityAnalysisResult, self).__init__(*args, **kwargs)
        self._is_essential = is_essential
        self._species_id = species_id
        self._exchange_fluxes = exchange_fluxes
        self._steps = steps
        self._biomass_fluxes = biomass_fluxes
        self._biomass = biomass
        self._variables = variables
        self._variables_fluxes = variables_fluxes

    @property
    def data_frame(self):

        data = [self._steps, self._exchange_fluxes]
        index = ["fraction", self._species_id]

        if self._biomass is not None:
            data.append(self._biomass_fluxes)
            index.append(self._biomass.id)

        if self._variables is not None and len(self._variables) > 0:
            for i, variable in enumerate(self._variables):
                data.append(self._variables_fluxes[:, i].tolist())
                if isinstance(variable, Reaction):
                    variable = variable.id
                index.append(variable)

        return DataFrame(data, index=index).T

    def plot(self, grid=None, width=None, height=None, *args, **kwargs):
        x_label = "Competition Level" if self._is_essential else 'Inhibition level'
        fig = figure(plot_width=width, plot_height=height, title=self._species_id,
                     x_axis_label=x_label, y_axis_label="Accumulation Level (mmol/gDW)",
                     toolbar_sticky=False)

        data = self.data_frame
        fig.extra_y_ranges = {}
        if self._biomass is not None:
            fig.extra_y_ranges["growth_rate"] = Range1d(start=0, end=data[self._biomass.id].max())
            fig.add_layout(LinearAxis(y_range_name="growth_rate", axis_label="Growth Rate [h^-1])"), 'right')
        if self._variables is not None and len(self._variables) > 0:
            fig.extra_y_ranges["flux"] = Range1d(start=min(data[data.columns[3:]].min()),
                                                 end=max(data[data.columns[3:]].max()))
            fig.add_layout(LinearAxis(y_range_name="flux", axis_label="Flux [mmol h^-1 gDW^-1])"), 'right')
        if self._biomass is None:
            fig.line(data['fraction'].apply(lambda v: v if self._is_essential else 1 - v) * 100,
                     data[self._species_id], line_color='orange')

        else:
            fig.line(data['fraction'].apply(lambda v: v if self._is_essential else 1 - v) * 100,
                     data[self._species_id], line_color='orange')
            fig.line(data['fraction'].apply(lambda v: v if self._is_essential else 1 - v) * 100,
                     data[self._biomass.id], line_color='green',
                     y_range_name="growth_rate")

        for var in self._variables:
            if isinstance(var, Reaction):
                var = var.id
            fig.line(data['fraction'].apply(lambda v: v if self._is_essential else 1 - v) * 100,
                     data[var], y_range_name="flux", legend=var)

        show(fig)


def sensitivity_analysis(model, metabolite, biomass=None, variables=None, is_essential=False, steps=10,
                         reference_dist=None, simulation_method=fba, **simulation_kwargs):

    if reference_dist is None:
        simulation_kwargs['reference'] = simulation_kwargs.get('reference', None) or pfba(model, objective=biomass)

    if variables is None:
        variables = []

    species_id = metabolite.id[:-2]

    exchange_fluxes = []
    biomass_fluxes = []
    variables_fluxes = []
    fractions = []

    for i, fraction in enumerate(frange(0, 1.1, steps)):
        variables_fluxes.append([])
        with TimeMachine() as tm:
            if is_essential:
                exchange = compete_metabolite(model, metabolite, simulation_kwargs['reference'], fraction,
                                              time_machine=tm)
            else:
                exchange = inhibit_metabolite(model, metabolite, simulation_kwargs['reference'], fraction,
                                              time_machine=tm)
            try:
                flux_dist = simulation_method(model, objective=biomass, **simulation_kwargs)

                if biomass is not None:
                    biomass_fluxes.append(flux_dist[biomass])

                for variable in variables:
                    variables_fluxes[i].append(flux_dist[variable])

                flux = flux_dist[exchange]

                exchange_fluxes.append(flux)
                fractions.append(fraction)
                logger.debug("Feasible: %s (%.3f) essential: %s flux: %.3f" %
                             (species_id, fraction, is_essential, flux))
            except Infeasible:
                logger.debug("Infeasible: %s (%.3f) essential: %s" % (species_id, fraction, is_essential))
                if biomass is not None:
                    biomass_fluxes.append(0)

                for _ in variables:
                    variables_fluxes[i].append(0)

                exchange_fluxes.append(0)
                fractions.append(fraction)

    return SensitivityAnalysisResult(species_id, exchange_fluxes, fractions, is_essential, biomass_fluxes, biomass,
                                     variables, array(variables_fluxes))
