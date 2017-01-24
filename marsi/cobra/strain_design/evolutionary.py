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
from IProgress import ProgressBar, Bar, Percentage
from cameo import fba, phenotypic_phase_plane, flux_variability_analysis
from cameo.core.strain_design import StrainDesignMethod
from cameo.exceptions import SolveError
from cameo.strain_design.heuristic.evolutionary.archives import ProductionStrainArchive
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_min_yield, \
    biomass_product_coupled_yield
from cameo.strain_design.heuristic.evolutionary_based import OptGeneResult
from cameo.util import TimeMachine
from cameo.visualization.plotting import plotter
from pandas import DataFrame
from sympy.tensor.tests.test_tensor import numpy

from marsi.cobra.metaheuristic.optimization import MetaboliteKnockoutOptimization

logger = logging.getLogger(__name__)

__all__ = ["OptGeneMet"]


class OptGeneMet(StrainDesignMethod):

    def __init__(self, model, evolutionary_algorithm=inspyred.ec.GA,
                 essential_metabolites=None, plot=True, *args, **kwargs):
        super(OptGeneMet, self).__init__(*args, **kwargs)
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
            inhibition_fraction=0.0, competition_fraction=0.5, simulation_method=fba, growth_coupled=False,
            max_evaluations=20000, population_size=200, time_machine=None, max_results=50, seed=None, **kwargs):
        """
        Parameters
        ----------
        target: str, Metabolite or Reaction
            The design target
        biomass: str, Metabolite or Reaction
            The biomass definition in the model
        substrate: str, Metabolite or Reaction
            The main carbon source
        max_knockouts: int
            Max number of knockouts allowed
        variable_size: bool
            If true, all candidates have the same size. Otherwise the candidate size can be from 1 to max_knockouts.
        inhibition_fraction: float, default: 0.0
            How much an antimetabolite inhibits the reactions.
        competition_fraction: float, default: 0.0
            How much an antimetabolite will increase the flux of the reactions producing the metabolite reactions.
        simulation_method: function
            Any method from cameo.flux_analysis.simulation or equivalent
        growth_coupled: bool
            If true will use the minimum flux rate to compute the fitness
        max_evaluations: int
            Number of evaluations before stop
        population_size: int
            Number of individuals in each generation
        time_machine: TimeMachine
            See TimeMachine
        max_results: int
            Max number of different designs to return if found.
        kwargs: dict
            Arguments for the simulation method.
        seed: int
            A seed for random.

        Returns
        -------
        OptGeneResult
        """

        target = self._model._reaction_for(target, time_machine=time_machine)
        biomass = self._model._reaction_for(biomass, time_machine=time_machine)
        substrate = self._model._reaction_for(substrate, time_machine=time_machine)

        if growth_coupled:
            objective_function = biomass_product_coupled_min_yield(biomass, target, substrate)
        else:
            objective_function = biomass_product_coupled_yield(biomass, target, substrate)

        optimization_algorithm = MetaboliteKnockoutOptimization(
            model=self._model,
            heuristic_method=self._algorithm,
            essential_metabolites=self._essential_metabolites,
            inhibition_fraction=inhibition_fraction,
            competition_fraction=competition_fraction,
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

        return OptGeneMetResult(inhibition_fraction, competition_fraction, optimization_algorithm.essential_metabolites,
                                self._model, result, objective_function, simulation_method, self._manipulation_type,
                                biomass, target, substrate, kwargs)


def process_metabolite_knockout_solution(model, solution, simulation_method, simulation_kwargs, inhibition_fraction,
                                         competition_fraction, essential_metabolites, biomass, target, substrate,
                                         objective_function, simplify=True):
    """

    Arguments
    ---------

    model: SolverBasedModel
        A constraint-based model
    solution: tuple - (metabolites, reactions)
        The output of a decoder
    simulation_method: function
        See see cameo.flux_analysis.simulation
    simulation_kwargs: dict
        Keyword arguments to run the simulation method
    inhibition_fraction: float
    competition_fraction: float
    essential_metabolites: float
    biomass: Reaction
        Cellular biomass reaction
    target: Reaction
        The strain design target
    substrate: Reaction
        The main carbon source uptake rate
    objective_function:
        A cameo.strain_design.heuristic.evolutionary.objective_functions.ObjectiveFunction
    cache: ProblemCache
        A problem cache for performance improvement

    Returns
    -------

    list
        A list with: metabolites, blocked_reactions, size, fva_min, fva_max, target flux, biomass flux, yield, fitness
    """

    fva_min, fva_max, target_flux, growth_rate, target_yield, *of = \
        _test_solution(model, solution, simulation_method, simulation_kwargs, inhibition_fraction,
                       competition_fraction, essential_metabolites, biomass, target, substrate, objective_function)

    tested = []
    original_solution = list(solution)
    while len(tested) < len(original_solution):
        to_remove = solution.pop(0)
        tested.append(to_remove)

        _fva_min, _fva_max, _target_flux, _growth_rate, _target_yield, *_of = \
            _test_solution(model, solution, simulation_method, simulation_kwargs, inhibition_fraction,
                           competition_fraction, essential_metabolites, biomass, target, substrate, objective_function)

        if _fva_min != fva_min or _growth_rate < growth_rate or _target_yield < target_yield or \
           any(_of[i] < of[i] for i in range(len(of))):
            solution.append(to_remove)

    if set(solution) != set(original_solution):
        fva_min, fva_max, target_flux, growth_rate, target_yield, *of = \
            _test_solution(model, solution, simulation_method, simulation_kwargs, inhibition_fraction,
                           competition_fraction, essential_metabolites, biomass, target, substrate, objective_function)

    return [solution, len(solution), fva_min, fva_max, target_flux, growth_rate, target_yield] + of


def _test_solution(model, solution, simulation_method, simulation_kwargs, inhibition_fraction, competition_fraction,
                   essential_metabolites, biomass, target, substrate, objective_function):
    """
    Returns
    -------

    fva_min, fva_max, target_flux, growth_rate, target_yield, of[0], [of[1], .., of[n-1]]
    """

    from marsi.processing.models import apply_anti_metabolite
    from marsi.utils import search_metabolites

    metabolites = [search_metabolites(model, species_id) for species_id in solution]
    with TimeMachine() as tm:
        for _metabolites in metabolites:
            apply_anti_metabolite(_metabolites, essential_metabolites, simulation_kwargs['reference'],
                                  inhibition_fraction=inhibition_fraction, competition_fraction=competition_fraction,
                                  time_machine=tm)

        reactions = objective_function.reactions
        flux_dist = simulation_method(model, reactions=reactions, objective=biomass, **simulation_kwargs)
        model.change_objective(biomass, time_machine=tm)

        fva = flux_variability_analysis(model, fraction_of_optimum=0.99, reactions=[target])
        target_yield = flux_dist[target] / abs(flux_dist[substrate])

        fitness = objective_function(model, flux_dist, solution)

        return fva.lower_bound(target), fva.upper_bound(target), flux_dist[target], flux_dist[biomass], target_yield, fitness


class OptGeneMetResult(OptGeneResult):
    __method_name__ = "OptGeneMet"

    def __init__(self, inhibition_fraction, competition_fraction, essential_metabolites, *args, **kwargs):
        super(OptGeneMetResult, self).__init__(*args, **kwargs)
        self.inhibition_fraction = inhibition_fraction
        self.competition_fraction = competition_fraction
        self.essential_metabolites = essential_metabolites

    @staticmethod
    def _generate_designs(knockouts, manipulation_type):
        return []

    @property
    def data_frame(self):
        if self._processed_solutions is None:
            self._process_solutions()

        return DataFrame(self._processed_solutions)

    def _process_solutions(self):
        processed_solutions = DataFrame(columns=["metabolites", "size", "fva_min", "fva_max",
                                                 "target_flux", "biomass_flux", "yield", "fitness"])

        if len(self._knockouts) == 0:
            logger.warn("No solutions found")
            self._processed_solutions = processed_solutions

        else:
            progress = ProgressBar(maxval=len(self._knockouts), widgets=["Processing solutions: ", Bar(), Percentage()])
            for i, solution in progress(enumerate(self._knockouts)):
                try:
                    processed_solutions.loc[i] = process_metabolite_knockout_solution(
                        self._model, solution, self._simulation_method, self._simulation_kwargs,
                        self.inhibition_fraction, self.competition_fraction, self.essential_metabolites,
                        self._biomass, self._target, self._substrate, self._objective_function)
                except SolveError as e:
                    logger.error(e)
                    processed_solutions.loc[i] = [numpy.nan for _ in processed_solutions.columns]

            self._processed_solutions = processed_solutions

    def display_on_map(self, index=0, map_name=None, palette="YlGnBu"):
        with TimeMachine() as tm:
            for ko in self.data_frame.loc[index, "metabolites"]:
                self._model.metabolites.get_by_id(ko).knock_out(tm)
            fluxes = self._simulation_method(self._model, **self._simulation_kwargs)
            fluxes.display_on_map(map_name=map_name, palette=palette)

    def plot(self, index=0, grid=None, width=None, height=None, title=None, palette=None, **kwargs):
        wt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])
        with TimeMachine() as tm:
            for ko in self.data_frame.loc[index, "metabolites"]:
                self._model.metabolites.get_by_id(ko).knock_out(tm)
            mt_production = phenotypic_phase_plane(self._model, objective=self._target, variables=[self._biomass])

        if title is None:
            title = "Production Envelope"

        dataframe = DataFrame(columns=["ub", "lb", "value", "strain"])
        for _, row in wt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "WT"]],
                            columns=dataframe.columns)
            dataframe = dataframe.append(_df)
        for _, row in mt_production.iterrows():
            _df = DataFrame([[row['objective_upper_bound'], row['objective_lower_bound'], row[self._biomass.id], "MT"]],
                            columns=dataframe.columns)
            dataframe = dataframe.append(_df)

        plot = plotter.production_envelope(dataframe, grid=grid, width=width, height=height, title=title,
                                           x_axis_label=self._biomass.id, y_axis_label=self._target.id, palette=palette)
        plotter.display(plot)
