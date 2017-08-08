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

from __future__ import absolute_import

import logging

from cameo.flux_analysis.analysis import find_essential_metabolites
from cameo.flux_analysis.simulation import pfba
from cameo.strain_design.heuristic.evolutionary.decoders import SetDecoder
from cameo.strain_design.heuristic.evolutionary.evaluators import TargetEvaluator
from cameo.strain_design.heuristic.evolutionary.optimization import TargetOptimization
from cobra.exceptions import OptimizationError

from marsi.cobra.flux_analysis.manipulation import knockout_metabolite
from marsi.utils import search_metabolites

logger = logging.getLogger(__name__)

__all__ = ["MetaboliteKnockoutOptimization"]

METABOLITE_KNOCKOUT_TYPE = "metabolite knockout"


class MetaboliteDecoder(SetDecoder):
    """
    Decoder for set representation. Converts an integer set into metabolites.

    Parameters
    ----------

    representation : list
        Reactions.
    model : SolverBasedModel

    """

    def __init__(self, representation, model, *args, **kwargs):
        super(MetaboliteDecoder, self).__init__(representation, model, *args, **kwargs)

    def __call__(self, individual, flat=False, **kwargs):
        """
        Parameters
        ----------

        **kwargs
        individual : list
            a list of integers
        flat : bool
            if True, returns strings. Otherwise returns Reaction.

        Returns
        -------
        list
            Decoded representation

        """
        metabolites = tuple(self.representation[index] for index in individual)
        return metabolites


class AntiMetaboliteEvaluator(TargetEvaluator):
    def __init__(self, essential_metabolites=None, inhibition_fraction=.0, competition_fraction=.0, *args, **kwargs):
        super(AntiMetaboliteEvaluator, self).__init__(*args, **kwargs)
        self.essential_metabolites = essential_metabolites or []
        self.inhibition_fraction = inhibition_fraction
        self.competition_fraction = competition_fraction

    def _evaluate_individual(self, individual):
        return self.evaluate_individual(individual)

    def evaluate_individual(self, individual):
        from marsi.cobra.flux_analysis.manipulation import apply_anti_metabolite

        specie_ids = self.decoder(individual)
        metabolite_targets = [search_metabolites(self.model, val) for val in specie_ids]
        with self.model:
            for metabolites in metabolite_targets:
                apply_anti_metabolite(self.model, metabolites, self.essential_metabolites,
                                      self.simulation_kwargs['reference'], self.inhibition_fraction,
                                      self.competition_fraction)

            try:
                solution = self.simulation_method(self.model, **self.simulation_kwargs)
                return self.objective_function(self.model, solution, metabolite_targets)
            except OptimizationError:
                return self.objective_function.worst_fitness()


class MetaboliteKnockoutEvaluator(TargetEvaluator):
    def __init__(self, essential_metabolites=None, *args, **kwargs):
        super(MetaboliteKnockoutEvaluator, self).__init__(*args, **kwargs)
        self.essential_metabolites = essential_metabolites or []

    def _evaluate_individual(self, individual):
        return self.evaluate_individual(individual)

    def evaluate_individual(self, individual):

        specie_ids = self.decoder(individual)
        metabolite_targets = [search_metabolites(self.model, val) for val in specie_ids]
        with self.model:
            for metabolites in metabolite_targets:
                for metabolite in metabolites:
                    knockout_metabolite(self.model, metabolite, ignore_transport=True, allow_accumulation=True)
            try:
                solution = self.simulation_method(self.model, **self.simulation_kwargs)
                return self.objective_function(self.model, solution, metabolite_targets)
            except OptimizationError:
                return self.objective_function.worst_fitness()


class MetaboliteKnockoutOptimization(TargetOptimization):
    """
    Knockout optimization using metabolites.

    Attributes
    ----------
    model : SolverBasedModel
        A constraint-based model.
    heuristic_method : inspyred.ec.EvolutionaryComputation
        An evolutionary algorithm.
    objective_function : objective function or list(objective function)
        The objectives for the algorithm to maximize.
    seed : int
        A seed for random. It is auto-generated if None is given.
    termination : inspyred.ec.terminators
        A termination criteria for the algorithm. The default is inspyred.ec.terminators.evaluation_termination.
    simulation_method: flux_analysis.simulation
        The method used to simulate the model.
    wt_reference: dict
        A reference initial state for the optimization. It is required for flux_analysis.simulation.lmoma and
        flux_analysis.simulation.room. If not given, it will be computed using flux_analysis.simulation.pfba
    metabolites: list
        A list of valid metabolites to knockout. If None, then all metabolites in the model will be knockout candidates
        except the ones defined in essential_metabolites
    essential_metabolites: list
        A list of metabolites that cannot be knocked out. If None, then all essential genes will be removed from the
        valid genes set.

    Methods
    -------
    run(view=config.default_view, maximize=True, **kwargs)


    See Also
    --------
    inspyred
    cameo.config.default_view

    Example
    -------
    >>> from cameo import models
    >>> model = models.bigg.iJO1366
    >>> from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
    >>> bpcy = biomass_product_coupled_yield(model.reactions.Ec_biomass_iJO1366_core_53p95,
    >>>                                      model.reactions.EX_succ_e,
    >>>                                      model.reactions.EX_glc__D_e)
    >>> knockout_optimization = MetaboliteKnockoutOptimization(model=model, objective_function=bpcy)
    >>> knockout_optimization.run(max_evaluations=50000)

    """
    def __init__(self, metabolites=None, essential_metabolites=None, n_carbons=2, compartments="c",
                 skip_essential_metabolites=False, *args, **kwargs):
        super(MetaboliteKnockoutOptimization, self).__init__(*args, **kwargs)

        if compartments is None:
            compartments = list(self.model.compartments.keys())
        self.compartments = compartments

        if metabolites is None:
            self.metabolites = set([met.id[:-2] for met in self.model.metabolites if
                                    met.elements.get('C', 0) >= n_carbons and met.compartment in self.compartments])
        else:
            self.metabolites = metabolites

        logger.debug("Computing essential reactions...")
        if skip_essential_metabolites:
            self.essential_metabolites = set()
        else:
            _essential_metabolites = find_essential_metabolites(self.model)
            self.essential_metabolites = set([m.id for m in _essential_metabolites])

        if essential_metabolites:
            self.essential_metabolites.update(essential_metabolites)

        self.representation = list(self.metabolites)
        self._target_type = METABOLITE_KNOCKOUT_TYPE
        self._decoder = MetaboliteDecoder(self.representation, self.model)
        self._evaluator = MetaboliteKnockoutEvaluator(
            model=self.model,
            decoder=self._decoder,
            objective_function=self.objective_function,
            simulation_method=self._simulation_method,
            simulation_kwargs=self._simulation_kwargs,
            essential_metabolites=self.essential_metabolites)

    @TargetOptimization.simulation_method.setter
    def simulation_method(self, simulation_method):
        if self._simulation_kwargs.get("reference", None) is None:
            logger.warning("No WT reference found, generating using pfba.")
            self._simulation_kwargs['reference'] = pfba(self.model).fluxes
            logger.warning("Reference successfully computed.")
        self._simulation_method = simulation_method

    @TargetOptimization.simulation_kwargs.setter
    def simulation_kwargs(self, simulation_kwargs):
        if self._simulation_method and simulation_kwargs.get("reference", None) is None:
            logger.warning("No WT reference found, generating using pfba.")
            simulation_kwargs['reference'] = pfba(self.model).fluxes
            logger.warning("Reference successfully computed.")
        self._simulation_kwargs = simulation_kwargs
