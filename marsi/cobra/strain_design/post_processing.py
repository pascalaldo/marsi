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

from IProgress import ProgressBar, Bar, Percentage
from cameo.core.target import ReactionKnockoutTarget, ReactionModulationTarget
from pandas import DataFrame

from cameo.util import TimeMachine
from cameo.exceptions import SolveError
from cameo.core.solver_based_model import SolverBasedModel
from cameo.core.reaction import Reaction
from cameo.strain_design.heuristic.evolutionary.objective_functions import ObjectiveFunction

from marsi.cobra.strain_design.target import AntiMetaboliteManipulationTarget, MetaboliteKnockoutTarget

logger = logging.getLogger(__name__)


class ReplaceKnockoutFinder(object):
    def __init__(self, model, simulation_method, simulation_kwargs, currency_metabolites):
        self.model = model
        self.simulation_method = simulation_method
        self.simulation_kwargs = simulation_kwargs
        self.currency_metabolites = currency_metabolites

    def __call__(self, targets, fitness, objective_function, max_loss=0.2):
        return test_reaction_knockout_replacements(self.model, targets, objective_function, fitness,
                                                   self.simulation_method, simulation_kwargs=self.simulation_kwargs,
                                                   ignore_metabolites=self.currency_metabolites, max_loss=max_loss)


def test_reaction_knockout_replacements(model, targets, objective_function, fitness_value, simulation_method,
                                        simulation_kwargs=None, reference=None, ignore_metabolites=None,
                                        ignore_transport=True, allow_accumulation=True, max_loss=0.2):

    if reference is None:
        reference = {}

    if simulation_kwargs is None:
        simulation_kwargs = {}

    assert isinstance(model, SolverBasedModel)
    assert isinstance(objective_function, ObjectiveFunction)

    def valid_loss(val, base):
        fitness_loss = fitness_value - val

        return fitness_loss / fitness_value < max_loss and val - base > 1e-6

    target_test_count = {test: 0 for test in targets}  # Keep track of which targets where tested
    possible_targets = list(targets)
    anti_metabolites = DataFrame(columns=['reaction', 'metabolites', 'base_fitness', 'fitness', 'loss'])
    index = 0

    def termination_criteria():
        logger.debug("Possible targets: %i" % len(possible_targets))
        logger.debug("Tested targets: %s" % target_test_count)
        logger.debug("Antimetabolites: %s" % str(anti_metabolites))
        return len(possible_targets) == 0 or all(count == 1 for count in target_test_count.values())

    # Stop when all targets have been replaced or tested more then once.
    while not termination_criteria():
        with TimeMachine() as tm:
            test = possible_targets.pop(0)
            target_test_count[test] += 1

            logger.debug("Testing target %s (%i times)" % (test, target_test_count[test]))

            reaction = model.reactions.get_by_id(test)
            for ko in possible_targets:
                model.reactions.get_by_id(ko).knock_out(tm)

            base_solution = simulation_method(model, **simulation_kwargs)
            base_fitness = objective_function(model, base_solution, targets)

            anti_metabolite_targets = find_anti_metabolite_knockouts(reaction,
                                                                     ref_flux=reference.get(test, 0),
                                                                     ignore_transport=ignore_transport,
                                                                     ignore_metabolites=ignore_metabolites,
                                                                     allow_accumulation=allow_accumulation)

            new_fitness_values = {}

            for species_id, target in anti_metabolite_targets.items():
                assert isinstance(target, AntiMetaboliteManipulationTarget)
                with TimeMachine() as another_tm:
                    target.apply(model, another_tm)
                    fitness = 0
                    try:
                        solution = simulation_method(model, **simulation_kwargs)
                        fitness = objective_function(model, solution, targets)
                    except SolveError:
                        logger.debug("Didn't solve %s" % species_id)
                        fitness = 0
                    finally:
                        if fitness not in new_fitness_values:
                            new_fitness_values[fitness] = []
                        new_fitness_values[fitness].append(species_id)

                        logger.debug("Applying %s yields %f (loss: %f)" %
                                     (species_id, fitness, (fitness_value - fitness)/fitness_value))

            new_fitness_values = {val: anti_mets for val, anti_mets in new_fitness_values.items()
                                  if valid_loss(val, base_fitness)}

            if len(new_fitness_values) == 0:
                # Put target back to test in other combinations.
                possible_targets.append(test)
                logger.debug("Return target %s, no replacement found" % test)
            else:
                max_fitness = max(new_fitness_values.keys())
                loss = fitness_value - max_fitness
                anti_metabolites.loc[index] = [test, new_fitness_values[max_fitness], base_fitness, max_fitness, loss]
                index += 1

    return anti_metabolites


def replace_knockouts(model, knockouts, objective_function, simulation_method, simulation_kwargs, currency_metabolites,
                      max_loss=0.2):
        knockout_replace_finder = ReplaceKnockoutFinder(model,
                                                        simulation_method=simulation_method,
                                                        simulation_kwargs=simulation_kwargs,
                                                        currency_metabolites=currency_metabolites)

        result = DataFrame(columns=['reactions', 'replaced_reaction', 'metabolites',
                                    'original_fitness', 'anti_metabolite_fitness'])
        pbar = ProgressBar(maxval=len(knockouts), widgets=["Replacing KOs: ", Bar(), Percentage()])
        replacements = [knockout_replace_finder(row.reactions, row.fitness, objective_function, max_loss)
                        for row in pbar(knockouts.itertuples(index=False))]
        i = 0
        for replace, original in zip(replacements, knockouts.itertuples(index=False)):
            if len(replace) == 0:
                continue

            for row in replace.itertuples(index=False):
                reactions = list(original.reactions)
                reactions.remove(row.reaction)
                for metabolite in row.metabolites:
                    result.loc[i] = [reactions, row.reaction, metabolite, original.fitness, row.fitness]
                    i += 1

        return result


def find_anti_metabolite_knockouts(reaction, ref_flux=0, ignore_metabolites=None, ignore_transport=True,
                                   allow_accumulation=True):
    """
    Generates a dictionary {species_id -> MetaboliteKnockoutTarget}.

    Parameters
    ----------
    reaction: cobra.Reaction
        A COBRA reaction
    ref_flux: float
        The flux from the reference state (0 if unkown)
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


def find_anti_metabolite_modulation(reaction, fold_change, ref_flux=0, ignore_metabolites=None,
                                    ignore_transport=True, allow_accumulation=True, essential_metabolites=None):
    """
    Generates a dictionary {species_id -> AntiMetaboliteManipulationTarget}.

    If the fold change > 1:
        1. Search for metabolites that are essential.
        2. Calculate fraction using a link function.

    If fold change < 1:
        1. Search for metabolites that are not essential
        2. Calculated fraction is 1 - fold change

    Parameters
    ----------
    reaction: cobra.Reaction
        A COBRA reaction.
    fold_change: float
        The fold change of the reaction flux.
    ref_flux: float
        The flux from the reference state (0 if unknown)
    ignore_metabolites: list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc.).
    ignore_transport: bool
        If False, also knockout the transport reactions.
    allow_accumulation: bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    essential_metabolites: list
        A list of essential metabolites.
    Returns
    -------
    dict



    Returns
    -------
    dict

    """

    assert isinstance(reaction, Reaction)
    assert isinstance(ignore_metabolites, (list, set, tuple))
    assert isinstance(essential_metabolites, (list, set, tuple))

    if fold_change > 1:
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
    if fold_change > 1:
        fraction = 10 * (1 - fold_change) / (1 + (1 - fold_change))
    else:
        fraction = 1 - fold_change

    for species_id in species_ids:
        result[species_id] = AntiMetaboliteManipulationTarget(species_id, fraction=fraction,
                                                              ignore_transport=ignore_transport,
                                                              allow_accumulation=allow_accumulation)

    return result


def convert_target(model, target, ignore_transport=True, ignore_metabolites=None,
                   allow_accumulation=True, reference=None):
    """
    Generates a dictionary {species_id -> MetaboliteKnockoutTarget}.

    Parameters
    ----------
    model: cobra.Model
        A COBRA model
    target: ReactionModulationTarget, ReactionKnockoutTarget
        The flux from the reference state (0 if unknown)
    ignore_metabolites: list
        A list of metabolites that should not be targeted (essential metabolites, currency metabolites, etc.)
    ignore_transport: bool
        If False, also knockout the transport reactions.
    allow_accumulation: bool
        If True, create an exchange reaction (unless already exists) to simulate accumulation of the metabolites.
    reference: dict
        A dictionary containing the flux values of a reference flux distribution.
    Returns
    -------
    dict
    """

    reference = reference or {}
    if isinstance(target, ReactionKnockoutTarget):

        substitutions = find_anti_metabolite_knockouts(target.get_model_target(model),
                                                       ref_flux=reference.get(target.id, 0),
                                                       ignore_transport=ignore_transport,
                                                       ignore_metabolites=ignore_metabolites,
                                                       allow_accumulation=allow_accumulation)

    elif isinstance(target, ReactionModulationTarget):
        substitutions = find_anti_metabolite_modulation(target.get_model_target(model),
                                                        target.fold_change,
                                                        ref_flux=reference.get(target.id, 0),
                                                        ignore_transport=ignore_transport,
                                                        ignore_metabolites=ignore_metabolites,
                                                        allow_accumulation=allow_accumulation)

    else:
        raise ValueError("Don't know what to do with the target of type %s" % type(target))

    return substitutions
