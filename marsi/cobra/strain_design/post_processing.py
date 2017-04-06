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

import logging

import numpy
from cameo.core.strain_design import StrainDesign
from cameo.core.target import ReactionKnockoutTarget, ReactionModulationTarget
from pandas import DataFrame

from cameo.util import TimeMachine
from cameo.exceptions import SolveError
from cameo.core.solver_based_model import SolverBasedModel
from cameo.core.reaction import Reaction
from cameo.strain_design.heuristic.evolutionary.objective_functions import ObjectiveFunction

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
        ignore_metabolites = []

    assert isinstance(reaction, Reaction)
    assert isinstance(ignore_metabolites, (list, set, tuple))
    assert isinstance(essential_metabolites, (list, set, tuple))

    if fold_change > 0:
        ignore_metabolites = list(ignore_metabolites) + essential_metabolites

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
        fraction = 1/(1 + numpy.exp(-0.5 * fold_change - 5))
    else:
        fraction = 1 - fold_change

    for species_id in species_ids:
        result[species_id] = AntiMetaboliteManipulationTarget(species_id, fraction=fraction,
                                                              ignore_transport=ignore_transport,
                                                              allow_accumulation=allow_accumulation)

    return result


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

    if simulation_kwargs is None:
        simulation_kwargs = {}

    replacements = DataFrame()
    for index, design in enumerate(results):
        with TimeMachine() as tm:
            design.apply(model, tm)
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
        essential_metabolites = []

    if ignore_metabolites is None:
        ignore_metabolites = []

    assert isinstance(model, SolverBasedModel)
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
    index = 0

    def termination_criteria():
        logger.debug("Targets: %i/%i" % (sum(target_test_count.values()), len(test_targets)))
        logger.debug("Anti metabolites: %s" % str(anti_metabolites))
        return len(test_targets) == 0 or all(count == 1 for count in target_test_count.values())

    # Stop when all targets have been replaced or tested more then once.
    while not termination_criteria():
        with TimeMachine() as tm:
            test_target = test_targets.pop(0)
            target_test_count[test_target.id] += 1

            logger.debug("Testing target %s" % test_target)
            assert test_target not in test_targets

            all_targets = test_targets + keep_targets

            for target in all_targets:
                target.apply(model, time_machine=tm)

            base_solution = simulation_method(model, **simulation_kwargs)
            base_fitness = objective_function(model, base_solution, test_targets)

            try:
                anti_metabolite_targets = convert_target(model, test_target, essential_metabolites,
                                                         ignore_transport=ignore_transport,
                                                         ignore_metabolites=ignore_metabolites,
                                                         allow_accumulation=allow_accumulation,
                                                         reference=reference)
                fitness2targets = {}

                for species_id, target in anti_metabolite_targets.items():
                    assert isinstance(target, AntiMetaboliteManipulationTarget)
                    with TimeMachine() as another_tm:
                        try:
                            target.apply(model, time_machine=another_tm, reference=reference)
                            new_solution = simulation_method(model, **simulation_kwargs)
                            new_fitness = objective_function(model, new_solution, all_targets)
                            logger.debug("New fitness %s" % new_fitness)
                            logger.debug("Solver objective value %s" % new_solution.objective_value)
                            for r in objective_function.reactions:
                                logger.debug("%s: %f" % (r, new_solution[r]))
                        except SolveError:
                            logger.debug("Cannot solve %s" % species_id)
                            new_fitness = 0
                        finally:
                            if new_fitness not in fitness2targets:
                                fitness2targets[new_fitness] = []
                            fitness2targets[new_fitness].append(target)
                            try:
                                logger.debug("Applying %s yields %f (loss: %f)" %
                                             (species_id, new_fitness, (new_fitness - fitness)/new_fitness))
                            except ZeroDivisionError:
                                logger.debug("Applying %s yields %f (loss: %f)" %
                                             (species_id, new_fitness, 1))

                # keep only targets that keep the fitness above the valid loss regarding the original fitness.
                fitness2targets = {fit: anti_mets for fit, anti_mets in fitness2targets.items()
                                   if valid_loss(fit, base_fitness)}

                if len(fitness2targets) == 0:
                    logger.debug("Return target %s, no replacement found" % test_target)
                else:
                    for fit, anti_mets in fitness2targets.items():
                        delta = fitness - fit
                        anti_metabolites.loc[index] = [StrainDesign(all_targets), test_target,
                                                       anti_mets, fitness, fit, delta]
                        index += 1
            except (ValueError, KeyError) as e:
                logger.error(str(e))
                continue
            finally:  # put the target back on the list.
                test_targets.append(test_target)

    return anti_metabolites


def convert_target(model, target, essential_metabolites, ignore_transport=True,
                   ignore_metabolites=None, allow_accumulation=True, reference=None):
    """
    Generates a dictionary {species_id -> MetaboliteKnockoutTarget}.

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

    Returns
    -------
    dict
    """

    if ignore_metabolites is None:
        ignore_metabolites = []

    if essential_metabolites is None:
        essential_metabolites = []

    reference = reference or {}
    if isinstance(target, ReactionKnockoutTarget):

        substitutions = find_anti_metabolite_knockouts(target.get_model_target(model),
                                                       ref_flux=reference.get(target.id, 0),
                                                       ignore_transport=ignore_transport,
                                                       ignore_metabolites=ignore_metabolites + essential_metabolites,
                                                       allow_accumulation=allow_accumulation)

    elif isinstance(target, ReactionModulationTarget):
        substitutions = find_anti_metabolite_modulation(target.get_model_target(model),
                                                        target.fold_change,
                                                        essential_metabolites,
                                                        ref_flux=reference.get(target.id, 0),
                                                        ignore_transport=ignore_transport,
                                                        ignore_metabolites=ignore_metabolites,
                                                        allow_accumulation=allow_accumulation)

    else:
        raise ValueError("Don't know what to do with the target of type %s" % type(target))

    return substitutions
