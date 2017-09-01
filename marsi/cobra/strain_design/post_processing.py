# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

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

import numpy
import six
from cameo.core.strain_design import StrainDesign
from cameo.core.target import ReactionKnockoutTarget, ReactionModulationTarget
from cameo.flux_analysis.structural import create_stoichiometric_array, find_coupled_reactions_nullspace
from cameo.flux_analysis.structural import nullspace
from cameo.strain_design.heuristic.evolutionary.objective_functions import ObjectiveFunction
from cobra.core.model import Model
from cobra.core.reaction import Reaction
from cobra.exceptions import OptimizationError
from optlang.exceptions import SolverError
from pandas import DataFrame, Series

from marsi.cobra import utils
from marsi.cobra.strain_design.target import AntiMetaboliteManipulationTarget, MetaboliteKnockoutTarget

logger = logging.getLogger(__name__)


def find_anti_metabolite_knockouts(reaction, ref_flux=0, ignore_metabolites=None, ignore_transport=True,
                                   allow_accumulation=True):
    """
    Generates a dictionary {species_id -> MetaboliteKnockoutTarget}.

    Parameters
    ----------
    reaction: cobra.Reaction
        A COBRA reaction
    ref_flux: float
        The flux from the reference state (0 if unknown)
    ignore_metabolites: list
        A list of metabolites that should not be targeted (currency metabolites, etc.)
    ignore_transport: bool
        If False, also knockout the transport reactions.
    allow_accumulation: bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.

    Returns
    -------
    dict

    """
    assert isinstance(reaction, Reaction)
    assert isinstance(ignore_metabolites, (list, set, tuple))

    if ref_flux != 0:
        substrates = [m for m, coefficient in reaction.metabolites.items()
                      if coefficient * ref_flux > 0 and m.id[:-2] not in ignore_metabolites]
    else:
        if reaction.reversibility:
            substrates = [m for m in reaction.metabolites.keys() if m.id[:-2] not in ignore_metabolites]
        else:
            substrates = [m for m, coefficient in reaction.metabolites.items()
                          if coefficient > 0 and m.id[:-2] not in ignore_metabolites]

    species_ids = [m.id[:-2] for m in substrates]
    result = {}
    for species_id in species_ids:
        result[species_id] = MetaboliteKnockoutTarget(species_id, ignore_transport, allow_accumulation)

    return result


def find_anti_metabolite_modulation(reaction, fold_change, essential_metabolites, ref_flux=0, ignore_metabolites=None,
                                    ignore_transport=True, allow_accumulation=True):
    """
    Generates a dictionary {species_id -> AntiMetaboliteManipulationTarget}.

    If the fold change > 0:
        1. Search for metabolites that are essential.
        2. Calculate fraction using a link function.

    If fold change < 0:
        1. Search for metabolites that are not essential
        2. Calculated fraction is 1 - fold change

    Parameters
    ----------
    reaction : cobra.Reaction
        A COBRA reaction.
    fold_change : float
        The fold change of the reaction flux.
    essential_metabolites : list
        A list of essential metabolites.
    ref_flux : float
        The flux from the reference state (0 if unknown)
    ignore_metabolites : list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc.).
    ignore_transport : bool
        If False, also knockout the transport reactions.
    allow_accumulation : bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.

    Returns
    -------
    dict
        {MetaboliteID -> AntiMetaboliteManipulationTarget}
    """

    if ignore_metabolites is None:
        ignore_metabolites = set(utils.CURRENCY_METABOLITES)

    if essential_metabolites is None and reaction.model is not None:
        essential_metabolites = utils.essential_species_ids(reaction.model)

    assert isinstance(reaction, Reaction)
    assert isinstance(ignore_metabolites, (list, set, tuple))
    assert isinstance(essential_metabolites, (list, set, tuple))

    if fold_change > 0:
        ignore_metabolites = set(ignore_metabolites) | set(essential_metabolites)

    if ref_flux != 0:
        substrates = [m for m, coefficient in reaction.metabolites.items()
                      if coefficient * ref_flux > 0 and m.id[:-2] not in ignore_metabolites]
    else:
        if reaction.reversibility:
            substrates = [m for m in reaction.metabolites.keys() if m.id[:-2] not in ignore_metabolites]
        else:
            substrates = [m for m, coefficient in reaction.metabolites.items()
                          if coefficient > 0 and m.id[:-2] not in ignore_metabolites]

    species_ids = [m.id[:-2] for m in substrates]
    result = {}

    # Use a link function to convert fold change into ]0, 1]
    if fold_change > 0:
        fraction = 1 / (1 + numpy.exp(-0.5 * fold_change - 5))
    else:
        fraction = 1 - fold_change

    for species_id in species_ids:
        result[species_id] = AntiMetaboliteManipulationTarget(species_id, fraction=fraction,
                                                              ignore_transport=ignore_transport,
                                                              allow_accumulation=allow_accumulation)

    return result


def convert_strain_design_results(model, results, objective_function, simulation_method, simulation_kwargs=None,
                                  ignore_metabolites=None, ignore_transport=True, allow_accumulation=True,
                                  nullspace_matrix=None, essential_metabolites=None, max_loss=0.2):
    """
    Converts a StrainDesignMethodResult into a DataFrame of possible substitutions.

    Parameters
    ----------
    model : cobra.Model
        A COBRA model.
    results : cameo.core.strain_design.StrainDesignMethodResult
        The results of a strain design method.
    objective_function : cameo.strain_design.heuristic.evolutionary.objective_functions.ObjectiveFunction
        The cellular objective to evaluate.
    simulation_method : cameo.flux_analysis.simulation.fba or equivalent
        The method to compute a flux distribution using a COBRA model.
    simulation_kwargs : dict
        The arguments for the simulation_method.
    ignore_metabolites : list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc).
    ignore_transport : bool
        If False, also knockout the transport reactions.
    allow_accumulation : bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    nullspace_matrix : numpy.ndarray
        The nullspace of the model.
    essential_metabolites : list
        A list of essential metabolites.
    max_loss : float
        A number between 0 and 1 for how much the fitness is allowed to drop with the metabolite target.

    Returns
    -------
    pandas.DataFrame
        A data frame with the possible replacements.
    """

    if essential_metabolites is None:
        essential_metabolites = utils.essential_species_ids(model)

    if ignore_metabolites is None:
        ignore_metabolites = set(utils.CURRENCY_METABOLITES)

    if nullspace_matrix is None:
        nullspace_matrix = nullspace(create_stoichiometric_array(model))

    if simulation_kwargs is None:
        simulation_kwargs = {}

    replacements = DataFrame()
    for index, design in enumerate(results):
        with model:
            design.apply(model)
            solution = simulation_method(model, **simulation_kwargs)
            fitness = objective_function(model, solution, design.targets)
        if fitness <= 0:
            continue

        res = convert_design(model, design, fitness, objective_function, simulation_method,
                             simulation_kwargs=simulation_kwargs, ignore_metabolites=ignore_metabolites,
                             ignore_transport=ignore_transport, allow_accumulation=allow_accumulation,
                             nullspace_matrix=nullspace_matrix, essential_metabolites=essential_metabolites,
                             max_loss=max_loss)

        res['index'] = index
        replacements = replacements.append(res, ignore_index=True)

    return replacements


def convert_design(model, strain_design, fitness, objective_function, simulation_method, simulation_kwargs=None,
                   ignore_metabolites=None, ignore_transport=True, allow_accumulation=True, nullspace_matrix=None,
                   essential_metabolites=None, max_loss=0.2):
    """
    Converts a StrainDesign into a DataFrame of possible substitutions.

    Parameters
    ----------
    model: cobra.Model
        A COBRA model.
    strain_design : cameo.core.strain_design.StrainDesign
        The results of a strain design method.
    objective_function : cameo.strain_design.heuristic.evolutionary.objective_functions.ObjectiveFunction
        The cellular objective to evaluate.
    simulation_method : cameo.flux_analysis.simulation.fba or equivalent
        The method to compute a flux distribution using a COBRA model.
    simulation_kwargs : dict
        The arguments for the simulation_method.
    ignore_metabolites : list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc).
    ignore_transport : bool
        If False, also knockout the transport reactions.
    allow_accumulation : bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    nullspace_matrix : numpy.ndarray
        The nullspace of the model.
    essential_metabolites : list
        A list of essential metabolites.
    max_loss : float
        A number between 0 and 1 for how much the fitness is allowed to drop with the metabolite target.

    Returns
    -------
    pandas.DataFrame
        A data frame with the possible replacements.
    """

    if simulation_kwargs is None:
        simulation_kwargs = {}

    if simulation_kwargs.get('reference', None) is None:
        reference = {}
    else:
        reference = simulation_kwargs['reference']

    if essential_metabolites is None:
        essential_metabolites = utils.essential_species_ids(model)

    if ignore_metabolites is None:
        ignore_metabolites = set(utils.CURRENCY_METABOLITES)

    if nullspace_matrix is None:
        nullspace_matrix = nullspace(create_stoichiometric_array(model))

    coupled_reactions = find_coupled_reactions_nullspace(model, ns=nullspace_matrix)
    coupled_reactions = [{r.id: c for r, c in six.iteritems(g)} for g in coupled_reactions]

    assert isinstance(model, Model)
    assert isinstance(objective_function, ObjectiveFunction)
    assert isinstance(nullspace_matrix, numpy.ndarray)

    def valid_loss(val, base):
        fitness_loss = fitness - val

        return fitness_loss / fitness < max_loss and val - base > 1e-6

    testable_targets = [t for t in strain_design.targets
                        if isinstance(t, (ReactionKnockoutTarget, ReactionModulationTarget))]

    non_testable_targets = {t for t in strain_design.targets
                            if not isinstance(t, (ReactionKnockoutTarget, ReactionModulationTarget))}
    coupled_targets = {}
    for t in testable_targets:
        group = next((g for g in coupled_reactions if t in g), None)
        if group:
            group = frozenset(group)
            if group not in coupled_targets:
                coupled_targets[group] = []
            coupled_targets[group].append(t)

    coupled_targets = [frozenset(t) for k, t in six.iteritems(coupled_targets) if len(t) > 1]
    non_coupled_targets = {t for t in testable_targets if not any(t in group for group in coupled_targets)}

    anti_metabolites = DataFrame(columns=['base_design', 'replaced_target', 'metabolite_targets',
                                          'old_fitness', 'fitness', 'delta'])

    logger.debug("Coupled groups: %i; non-coupled targets %i" % (len(coupled_targets), len(non_coupled_targets)))

    with model as base_model:
        for target in non_testable_targets:
            target.apply(model)

        # test non-coupled targets as before
        for target in non_coupled_targets:
            with base_model:
                for other_target in non_coupled_targets:
                    if other_target != target:
                        other_target.apply(base_model)

                anti_metabolite_targets = convert_target(base_model, target, essential_metabolites,
                                                         ignore_metabolites=ignore_metabolites,
                                                         ignore_transport=ignore_transport,
                                                         allow_accumulation=allow_accumulation,
                                                         reference=reference)

                base_solution = simulation_method(base_model, **simulation_kwargs)
                base_fitness = objective_function(base_model, base_solution, strain_design.targets)

                test_target_substitutions(base_model, strain_design.targets, target, anti_metabolite_targets,
                                          objective_function, fitness, base_fitness, simulation_method,
                                          simulation_kwargs, reference, valid_loss, anti_metabolites)

        # test coupled targets as whole
        for target_group in coupled_targets:
            with base_model:

                for other_target_group in non_coupled_targets:
                    if other_target_group != target_group:
                        for other_target in other_target_group:
                            other_target.apply(model)

                anti_metabolite_targets = convert_target_group(base_model, target_group, essential_metabolites,
                                                               ignore_transport=ignore_transport,
                                                               allow_accumulation=allow_accumulation,
                                                               reference=reference)

                base_solution = simulation_method(model, **simulation_kwargs)
                base_fitness = objective_function(model, base_solution, strain_design.targets)

                test_target_substitutions(base_model, strain_design.targets, target_group, anti_metabolite_targets,
                                          objective_function, fitness, base_fitness, simulation_method,
                                          simulation_kwargs, reference, valid_loss, anti_metabolites)

        return anti_metabolites


def test_target_substitutions(model, all_targets, target, replacement_targets, objective_function, fitness,
                              base_fitness, simulation_method, simulation_kwargs, reference, loss_validation, results):
    fitness2targets = {}
    index = len(results)
    for species_id, replacement_target in replacement_targets.items():
        assert isinstance(replacement_target, AntiMetaboliteManipulationTarget)
        with model:
            try:
                replacement_target.apply(model, reference=reference)
                new_solution = simulation_method(model, **simulation_kwargs)
                new_fitness = objective_function(model, new_solution, all_targets)
                logger.debug("New fitness %s" % new_fitness)
                logger.debug("Solver objective value %s" % new_solution.objective_value)
                for r in objective_function.reactions:
                    logger.debug("%s: %f" % (r, new_solution[r]))
            except OptimizationError:
                logger.debug("Cannot solve %s" % species_id)
                new_fitness = 0
            finally:
                if new_fitness not in fitness2targets:
                    fitness2targets[new_fitness] = []
                fitness2targets[new_fitness].append(replacement_target)
                try:
                    logger.debug("Applying %s yields %f (loss: %f)" %
                                 (species_id, new_fitness, (new_fitness - fitness) / new_fitness))
                except ZeroDivisionError:
                    logger.debug("Applying %s yields %f (loss: %f)" % (species_id, new_fitness, 1))

        # keep only targets that keep the fitness above the valid loss regarding the original fitness.
        fitness2targets = {fit: anti_mets for fit, anti_mets in fitness2targets.items()
                           if loss_validation(fit, base_fitness)}

        if len(fitness2targets) == 0:
            logger.debug("Return target %s, no replacement found" % target)
        else:
            for fit, anti_mets in fitness2targets.items():
                delta = fitness - fit
                results.loc[index] = [StrainDesign(all_targets), target, tuple(anti_mets), fitness, fit, delta]
                index += 1


def replace_strain_design_results(model, results, objective_function, simulation_method, simulation_kwargs=None,
                                  ignore_metabolites=None, ignore_transport=True, allow_accumulation=True,
                                  essential_metabolites=None, max_loss=0.2):
    """
    Converts a StrainDesignMethodResult into a DataFrame of possible substitutions.

    Parameters
    ----------
    model : cobra.Model
        A COBRA model.
    results : cameo.core.strain_design.StrainDesignMethodResult
        The results of a strain design method.
    objective_function : cameo.strain_design.heuristic.evolutionary.objective_functions.ObjectiveFunction
        The cellular objective to evaluate.
    simulation_method : cameo.flux_analysis.simulation.fba or equivalent
        The method to compute a flux distribution using a COBRA model.
    simulation_kwargs : dict
        The arguments for the simulation_method.
    ignore_metabolites : list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc).
    ignore_transport : bool
        If False, also knockout the transport reactions.
    allow_accumulation : bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    essential_metabolites : list
        A list of essential metabolites.
    max_loss : float
        A number between 0 and 1 for how much the fitness is allowed to drop with the metabolite target.

    Returns
    -------
    pandas.DataFrame
        A data frame with the possible replacements.
    """

    if essential_metabolites is None:
        essential_metabolites = utils.essential_species_ids(model)

    if ignore_metabolites is None:
        ignore_metabolites = set(utils.CURRENCY_METABOLITES)

    if simulation_kwargs is None:
        simulation_kwargs = {}

    replacements = DataFrame()
    for index, design in enumerate(results):
        with model:
            design.apply(model)
            solution = simulation_method(model, **simulation_kwargs)
            fitness = objective_function(model, solution, design.targets)
        if fitness <= 0:
            continue

        res = replace_design(model, design, fitness, objective_function, simulation_method,
                             simulation_kwargs=simulation_kwargs, ignore_metabolites=ignore_metabolites,
                             ignore_transport=ignore_transport, allow_accumulation=allow_accumulation,
                             essential_metabolites=essential_metabolites, max_loss=max_loss)

        res['index'] = index
        replacements = replacements.append(res, ignore_index=True)

    return replacements


def replace_design(model, strain_design, fitness, objective_function, simulation_method, simulation_kwargs=None,
                   ignore_metabolites=None, ignore_transport=True, allow_accumulation=True,
                   essential_metabolites=None, max_loss=0.2, allow_modulation=True):
    """
    Converts a StrainDesign into a DataFrame of possible substitutions.

    Parameters
    ----------
    model: cobra.Model
        A COBRA model.
    strain_design : cameo.core.strain_design.StrainDesign
        The results of a strain design method.
    objective_function : cameo.strain_design.heuristic.evolutionary.objective_functions.ObjectiveFunction
        The cellular objective to evaluate.
    simulation_method : cameo.flux_analysis.simulation.fba or equivalent
        The method to compute a flux distribution using a COBRA model.
    simulation_kwargs : dict
        The arguments for the simulation_method.
    ignore_metabolites : list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc).
    ignore_transport : bool
        If False, also knockout the transport reactions.
    allow_accumulation : bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    essential_metabolites : list
        A list of essential metabolites.
    max_loss : float
        A number between 0 and 1 for how much the fitness is allowed to drop with the metabolite target.
    allow_modulation : bool
        If False does not allow modulation targets (MILP and QP are not compatible)

    Returns
    -------
    pandas.DataFrame
        A data frame with the possible replacements.
    """

    if simulation_kwargs is None:
        simulation_kwargs = {}

    if simulation_kwargs.get('reference', None) is None:
        reference = {}
    else:
        reference = simulation_kwargs['reference']

    if essential_metabolites is None:
        essential_metabolites = utils.essential_species_ids(model)

    if ignore_metabolites is None:
        ignore_metabolites = set(utils.CURRENCY_METABOLITES)

    assert isinstance(model, Model)
    assert isinstance(objective_function, ObjectiveFunction)

    def valid_loss(val, base):
        fitness_loss = fitness - val

        return fitness_loss / fitness < max_loss and val - base > 1e-6

    # Keep track of which targets where tested
    target_test_count = {test.id: 0 for test in strain_design.targets if isinstance(test, ReactionModulationTarget)}
    test_targets = [t for t in strain_design.targets if isinstance(t, ReactionModulationTarget)]
    keep_targets = [t for t in strain_design.targets if not isinstance(t, ReactionModulationTarget)]
    anti_metabolites = DataFrame(columns=['base_design', 'replaced_target', 'metabolite_targets',
                                          'old_fitness', 'fitness', 'delta'])

    def termination_criteria():
        logger.debug("Targets: %i/%i" % (sum(target_test_count.values()), len(test_targets)))
        logger.debug("Anti metabolites: %s" % str(anti_metabolites))
        return len(test_targets) == 0 or all(count == 1 for count in target_test_count.values())

    # Stop when all targets have been replaced or tested more then once.
    while not termination_criteria():
        with model as base_model:
            test_target = test_targets.pop(0)
            target_test_count[test_target.id] += 1

            logger.debug("Testing target %s" % test_target)
            assert test_target not in test_targets

            all_targets = test_targets + keep_targets

            for target in all_targets:
                target.apply(model)

            base_solution = simulation_method(base_model, **simulation_kwargs)
            base_fitness = objective_function(base_model, base_solution, test_targets)

            try:
                anti_metabolite_targets = convert_target(base_model, test_target, essential_metabolites,
                                                         ignore_transport=ignore_transport,
                                                         ignore_metabolites=ignore_metabolites,
                                                         allow_accumulation=allow_accumulation,
                                                         reference=reference,
                                                         allow_modulation=allow_modulation)

                test_target_substitutions(base_model, all_targets, test_target, anti_metabolite_targets,
                                          objective_function, fitness, base_fitness, simulation_method,
                                          simulation_kwargs, reference, valid_loss, anti_metabolites)
            except (ValueError, KeyError, SolverError) as e:
                logger.error(str(e))
                continue
            finally:  # put the target back on the list.
                test_targets.append(test_target)

    anti_metabolites.drop_duplicates(['replaced_target', 'metabolite_targets'], inplace=True)
    anti_metabolites.index = [i for i in range(len(anti_metabolites))]

    return anti_metabolites


def convert_target(model, target, essential_metabolites, ignore_transport=True, ignore_metabolites=None,
                   allow_accumulation=True, reference=None, allow_modulation=True):
    """
    Generates a dictionary {species_id -> (MetaboliteKnockoutTarget, MetaboliteModulationTarget)}.

    Parameters
    ----------
    model : cobra.Model
        A COBRA model
    target : ReactionModulationTarget, ReactionKnockoutTarget
        The flux from the reference state (0 if unknown)
    ignore_metabolites : list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc.)
    ignore_transport : bool
        If False, also knockout the transport reactions.
    allow_accumulation : bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    reference : dict
        A dictionary containing the flux values of a reference flux distribution.
    essential_metabolites : list
        A list of essential metabolites
    allow_modulation : bool
        If False does not allow modulation targets (MILP and QP are not compatible).

    Returns
    -------
    dict
    """

    if ignore_metabolites is None:
        ignore_metabolites = set(utils.CURRENCY_METABOLITES)

    if essential_metabolites is None:
        essential_metabolites = utils.essential_species_ids(model)

    if reference is None:
        reference = {}
    elif isinstance(reference, Series):
        reference = reference.to_dict()

    if isinstance(target, ReactionKnockoutTarget):

        substitutions = find_anti_metabolite_knockouts(target.get_model_target(model),
                                                       ref_flux=reference.get(target.id, 0),
                                                       ignore_transport=ignore_transport,
                                                       ignore_metabolites=ignore_metabolites | essential_metabolites,
                                                       allow_accumulation=allow_accumulation)

    elif isinstance(target, ReactionModulationTarget) and allow_modulation:
        substitutions = find_anti_metabolite_modulation(target.get_model_target(model),
                                                        target.fold_change,
                                                        essential_metabolites,
                                                        ref_flux=reference.get(target.id, 0),
                                                        ignore_transport=ignore_transport,
                                                        ignore_metabolites=ignore_metabolites,
                                                        allow_accumulation=allow_accumulation)

    else:
        substitutions = {}

    return substitutions


def convert_target_group(model, target_group, essential_metabolites, ignore_transport=True, ignore_metabolites=None,
                         allow_accumulation=True, reference=None, allow_modulation=True):
    """
    Generates a dictionary {species_id -> MetaboliteKnockoutTarget}.

    Parameters
    ----------
    model : cobra.Model
        A COBRA model
    target_group : dict
        A dictionary of (target: relative coefficient)
    ignore_metabolites : set
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc.)
    ignore_transport : bool
        If False, also knockout the transport reactions.
    allow_accumulation : bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    reference : dict
        A dictionary containing the flux values of a reference flux distribution.
    essential_metabolites : set
        A list of essential metabolites.
    allow_modulation : bool
        If False does not allow modulation targets (MILP and QP are not compatible).

    Returns
    -------
    tuple
        A dictionary with the antimetabolites to test and another dictionary with the reaction where those metabolites
        belong.
    """

    if ignore_metabolites is None:
        ignore_metabolites = set(utils.CURRENCY_METABOLITES)

    if essential_metabolites is None:
        essential_metabolites = utils.essential_species_ids(model)

    if reference is None:
        reference = {}
    elif isinstance(reference, Series):
        reference = reference.to_dict()

    substitutions = {}
    metabolites_targets = {}
    for target in target_group:
        if isinstance(target, ReactionKnockoutTarget):
            ignore = ignore_metabolites | essential_metabolites
            _substitutions = find_anti_metabolite_knockouts(target.get_model_target(model),
                                                            ref_flux=reference.get(target.id, 0),
                                                            ignore_transport=ignore_transport,
                                                            ignore_metabolites=ignore,
                                                            allow_accumulation=allow_accumulation)

        elif isinstance(target, ReactionModulationTarget) and allow_modulation:
            _substitutions = find_anti_metabolite_modulation(target.get_model_target(model),
                                                             target.fold_change,
                                                             essential_metabolites,
                                                             ref_flux=reference.get(target.id, 0),
                                                             ignore_transport=ignore_transport,
                                                             ignore_metabolites=ignore_metabolites,
                                                             allow_accumulation=allow_accumulation)
        else:
            _substitutions = {}

        substitutions.update(_substitutions)
        for substitution in _substitutions:
            if substitution not in metabolites_targets:
                metabolites_targets[substitution] = set()
            metabolites_targets[substitution].append(target)

    return substitutions, metabolites_targets
