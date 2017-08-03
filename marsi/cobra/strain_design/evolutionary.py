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

import inspyred
import numpy
from IProgress import ProgressBar, Bar, Percentage
from cameo.core.strain_design import StrainDesignMethod, StrainDesign, StrainDesignMethodResult
from cameo.core.utils import get_reaction_for
from cameo.flux_analysis.analysis import phenotypic_phase_plane, flux_variability_analysis
from cameo.flux_analysis.simulation import fba
from cameo.strain_design.heuristic.evolutionary.archives import ProductionStrainArchive
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_min_yield, \
    biomass_product_coupled_yield
from cameo.visualization.plotting import plotter
from cobra.core import Model
from cobra.exceptions import OptimizationError
from pandas import DataFrame

from marsi.cobra.flux_analysis.manipulation import knockout_metabolite
from marsi.cobra.strain_design.metaheuristic import MetaboliteKnockoutOptimization
from marsi.cobra.strain_design.target import MetaboliteKnockoutTarget
from marsi.utils import search_metabolites

logger = logging.getLogger(__name__)

__all__ = ["OptMet"]


class OptMet(StrainDesignMethod):

    def __init__(self, model, evolutionary_algorithm=inspyred.ec.GA,
                 essential_metabolites=None, plot=True, *args, **kwargs):
        super(OptMet, self).__init__(*args, **kwargs)
        self._model = model
        self._algorithm = evolutionary_algorithm
        self._optimization_algorithm = None
        self._essential_metabolites = essential_metabolites
        self._plot = plot
        self._manipulation_type = "metabolites"

    @property
    def manipulation_type(self):
        return self._manipulation_type

    @property
    def plot(self):
        return self._plot

    @plot.setter
    def plot(self, plot):
        self._plot = plot
        if self._optimization_algorithm is not None:
            self._optimization_algorithm.plot = plot

    def run(self, target=None, biomass=None, substrate=None, max_knockouts=5, variable_size=True,
            simulation_method=fba, growth_coupled=False, max_evaluations=20000, population_size=200,
            max_results=50, seed=None, **kwargs):
        """
        Parameters
        ----------
        target : str, Metabolite or Reaction
            The design target
        biomass : str, Metabolite or Reaction
            The biomass definition in the model
        substrate : str, Metabolite or Reaction
            The main carbon source
        max_knockouts : int
            Max number of knockouts allowed
        variable_size : bool
            If true, all candidates have the same size. Otherwise the candidate size can be from 1 to max_knockouts.
        simulation_method: function
            Any method from cameo.flux_analysis.simulation or equivalent
        growth_coupled : bool
            If true will use the minimum flux rate to compute the fitness
        max_evaluations : int
            Number of evaluations before stop
        population_size : int
            Number of individuals in each generation
        max_results : int
            Max number of different designs to return if found.
        kwargs : dict
            Arguments for the simulation method.
        seed : int
            A seed for random.

        Returns
        -------
        OptMetResult
        """

        target = get_reaction_for(self._model, target)
        biomass = get_reaction_for(self._model, biomass)
        substrate = get_reaction_for(self._model, substrate)

        if growth_coupled:
            objective_function = biomass_product_coupled_min_yield(biomass, target, substrate)
        else:
            objective_function = biomass_product_coupled_yield(biomass, target, substrate)

        optimization_algorithm = MetaboliteKnockoutOptimization(
            model=self._model,
            heuristic_method=self._algorithm,
            essential_metabolites=self._essential_metabolites,
            objective_function=objective_function,
            plot=self.plot)

        optimization_algorithm.simulation_kwargs = kwargs
        optimization_algorithm.simulation_method = simulation_method
        optimization_algorithm.archiver = ProductionStrainArchive()

        result = optimization_algorithm.run(max_evaluations=max_evaluations,
                                            pop_size=population_size,
                                            max_size=max_knockouts,
                                            variable_size=variable_size,
                                            maximize=True,
                                            max_archive_size=max_results,
                                            seed=seed,
                                            **kwargs)

        kwargs.update(optimization_algorithm.simulation_kwargs)
        return OptMetResult(self._model, result, objective_function, simulation_method, self._manipulation_type,
                            biomass, target, substrate, kwargs)


class OptMetResult(StrainDesignMethodResult):
    __method_name__ = "OptMet"

    def __init__(self, model, knockouts, objective_function, simulation_method, manipulation_type, biomass, target,
                 substrate, simulation_kwargs, *args, **kwargs):

        assert isinstance(model, Model)

        self._model = model
        self._knockouts = knockouts
        self._objective_function = objective_function
        self._simulation_method = simulation_method
        self._manipulation_type = manipulation_type
        self._biomass = biomass
        self._target = target
        self._substrate = substrate
        self._processed_solutions = None
        self._simulation_kwargs = simulation_kwargs

        super(OptMetResult, self).__init__(self.designs, *args, **kwargs)

    @property
    def designs(self):
        if self._processed_solutions is None:
            self._process_solutions()
        return self._processed_solutions.designs.tolist()

    def _repr_html_(self):
        return """
        <h3>OptMet Result</h3>
        <ul>
            <li>Simulation: %s<br/></li>
            <li>Objective Function: %s<br/></li>
        </ul>
        %s
        """ % (self._simulation_method.__name__,
               self._objective_function._repr_latex_(),
               self.data_frame._repr_html_())

    @property
    def data_frame(self):
        if self._processed_solutions is None:
            self._process_solutions()

        return DataFrame(self._processed_solutions)

    def _process_solutions(self):
        processed_solutions = DataFrame(columns=["designs", "size", "fva_min", "fva_max",
                                                 "target_flux", "biomass_flux", "yield", "fitness"])

        if len(self._knockouts) == 0:
            logger.warning("No solutions found")
            self._processed_solutions = processed_solutions

        else:
            progress = ProgressBar(maxval=len(self._knockouts), widgets=["Processing solutions: ", Bar(), Percentage()])
            for i, solution in progress(enumerate(self._knockouts)):
                try:
                    processed_solutions.loc[i] = process_metabolite_knockout_solution(
                        self._model, solution, self._simulation_method, self._simulation_kwargs,
                        self._biomass, self._target, self._substrate, self._objective_function)
                except OptimizationError as e:
                    logger.error(e)
                    processed_solutions.loc[i] = [numpy.nan for _ in processed_solutions.columns]

            self._processed_solutions = processed_solutions

    def display_on_map(self, index=0, map_name=None, palette="YlGnBu"):
        with self._model:
            for ko in self.data_frame.loc[index, "metabolites"]:
                knockout_metabolite(self._model, self._model.metabolites.get_by_id(ko))
            fluxes = self._simulation_method(self._model, **self._simulation_kwargs)
            fluxes.display_on_map(map_name=map_name, palette=palette)

    def plot(self, index=0, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        wt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])
        with self._model:
            for ko in self.data_frame.loc[index, "metabolites"]:
                knockout_metabolite(self._model, self._model.metabolites.get_by_id(ko))
            mt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])

        if title is None:
            title = "Production Envelope"

        data_frame = DataFrame(columns=["ub", "lb", "value", "strain"])
        for _, row in wt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "WT"]],
                            columns=data_frame.columns)
            data_frame = data_frame.append(_df)
        for _, row in mt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "MT"]],
                            columns=data_frame.columns)
            data_frame = data_frame.append(_df)

        plot = plotter.production_envelope(data_frame, grid=grid, width=width, height=height, title=title,
                                           x_axis_label=self._biomass.id, y_axis_label=self._target.id, palette=palette)
        plotter.display(plot)


def process_metabolite_knockout_solution(model, solution, simulation_method, simulation_kwargs, biomass, target,
                                         substrate, objective_function):
    """
    Parameters
    ----------

    model : SolverBasedModel
        A constraint-based model
    solution : tuple
        The output of a decoder (metabolite species)
    simulation_method : function
        See see cameo.flux_analysis.simulation
    simulation_kwargs : dict
        Keyword arguments to run the simulation method
    biomass : Reaction
        Cellular biomass reaction
    target : Reaction
        The strain design target
    substrate : Reaction
        The main carbon source uptake rate
    objective_function : function
        A cameo.strain_design.heuristic.evolutionary.objective_functions.ObjectiveFunction

    Returns
    -------

    list
        A list with: metabolites, size, fva_min, fva_max, target flux, biomass flux, yield, fitness
    """

    metabolite_targets = [search_metabolites(model, species_id) for species_id in solution]
    with model:
        for metabolites in metabolite_targets:
            for metabolite in metabolites:
                knockout_metabolite(model, metabolite)

        reactions = objective_function.reactions
        flux_dist = simulation_method(model, reactions=reactions, objective=biomass, **simulation_kwargs)
        model.objective = biomass

        fva = flux_variability_analysis(model, fraction_of_optimum=0.99, reactions=[target])
        target_yield = flux_dist[target] / abs(flux_dist[substrate])

        fitness = objective_function(model, flux_dist, solution)

        fva_min, fva_max = fva.lower_bound(target), fva.upper_bound(target)
        target_flux, biomass_flux = flux_dist[target], flux_dist[biomass]

        design = StrainDesign(MetaboliteKnockoutTarget(species_id=species_id) for species_id in solution)

        return design, len(design), fva_min, fva_max, target_flux, biomass_flux, target_yield, fitness
