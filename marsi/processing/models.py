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

from numbers import Number

from cameo.core import SolverBasedModel
from cameo.exceptions import SolveError
from cameo.flux_analysis.simulation import FluxDistributionResult
from cameo.strain_design.heuristic.evolutionary.objective_functions import ObjectiveFunction
from cameo.util import TimeMachine
from cobra import Reaction
from functools import lru_cache

from collections import Counter
from pandas import DataFrame

from marsi import bigg_api
from marsi.io.bigg import bigg_metabolites
from marsi.mining.enrichment import inchi_from_chebi, inchi_from_kegg

__all__ = ['find_inchi_for_bigg_metabolite', 'inhibit_metabolite', 'compete_metabolite',
           'degree', 'in_degree', 'out_degree']

logger = logging.getLogger(__name__)

DATABASE_LINKS = 'database_links'

CHEBI = 'CHEBI'
KEGG = "KEGG Compound"


@lru_cache(maxsize=1024)
def find_inchi_for_bigg_metabolite(model, metabolite_id):
    try:
        links = bigg_metabolites.loc[metabolite_id].database_links
    except KeyError:
        metabolite_data = bigg_api.get_model_metabolite(model.id, metabolite_id)
        links = metabolite_data[DATABASE_LINKS]
    inchi_keys = []
    if CHEBI in links:
        inchi_keys += [inchi_from_chebi(link['id']) for link in links[CHEBI]]

    if KEGG in links:
        inchi_keys += [inchi_from_kegg(link['id']) for link in links[KEGG]]

    counter = Counter(inchi_keys)
    del counter[None]

    if len(counter) > 0:
        return max(counter, key=lambda key: counter[key])
    else:
        raise ValueError(metabolite_id)


def degree(metabolite):
    """
    Degree of a metabolite is the number of reactions coming in and out.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A Metabolite.

    Returns
    -------
    int
        The degree.

    """
    return sum(2 if r.reversibility else 1 for r in metabolite.reactions)


def in_degree(metabolite):
    """
    In degree of a metabolite is the number of reactions coming in.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A Metabolite.

    Returns
    -------
    int
        The degree.

    """
    return len([r for r in metabolite.reactions if r.metabolites[metabolite] < 0 or r.reversibility])


def out_degree(metabolite):
    """
    Degree of a metabolite is the number of reactions coming out.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A Metabolite.

    Returns
    -------
    int
        The degree.

    """
    return len([r for r in metabolite.reactions if r.metabolites[metabolite] > 0 or r.reversibility])


def apply_antimetabolite(metabolites, essential_metabolites, reference, inhibition_fraction=.0, competition_fraction=.0,
                         allow_accumulation=True, ignore_transport=True, time_machine=None):
    exchanges = set()
    if any(met.id in essential_metabolites for met in metabolites):
        for metabolite in metabolites:
            exchanges.add(compete_metabolite(metabolite,
                                             reference,
                                             allow_accumulation=allow_accumulation,
                                             ignore_transport=ignore_transport,
                                             fraction=competition_fraction,
                                             time_machine=time_machine))
    else:
        for metabolite in metabolites:
            exchanges.add(inhibit_metabolite(metabolite,
                                             reference,
                                             allow_accumulation=allow_accumulation,
                                             ignore_transport=ignore_transport,
                                             fraction=inhibition_fraction,
                                             time_machine=time_machine))

    return exchanges


def search_metabolites(model, species_id):
    return model.metabolites.query(lambda mid: mid[:-2] == species_id, attribute='id')


def test_reaction_knockout_replacements(model, targets, objective_function, fitness_value, simulation_method,
                                        simulation_kwargs=None, reference=None, ignore_metabolites=None,
                                        ignore_transport=True, allow_accumulation=True, max_lost=0.2):

    if reference is None:
        reference = {}

    if simulation_kwargs is None:
        simulation_kwargs = {}

    assert isinstance(model, SolverBasedModel)
    assert isinstance(objective_function, ObjectiveFunction)

    def valid_loss(val, base):
        fitness_loss = fitness_value - val

        return fitness_loss / fitness_value < max_lost and val - base > 1e-6

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

            anti_metabolite_targets = knockout_antimetabolite_substitution(model,
                                                                           reaction,
                                                                           ref_flux=reference.get(test, 0),
                                                                           ignore_transport=ignore_transport,
                                                                           ignore_metabolites=ignore_metabolites,
                                                                           allow_accumulation=True)

            new_fitness_values = {}

            for anti_metabolite, apply in anti_metabolite_targets.items():
                with TimeMachine() as another_tm:
                    apply(another_tm)
                    fitness = 0
                    try:
                        solution = simulation_method(model, **simulation_kwargs)
                        fitness = objective_function(model, solution, targets)
                    except SolveError:
                        logger.debug("Didn't solve %s" % anti_metabolite)
                        fitness = 0
                    finally:
                        if fitness not in new_fitness_values:
                            new_fitness_values[fitness] = []
                        new_fitness_values[fitness].append(anti_metabolite)

                        logger.debug("Applying %s yields %f (loss: %f)" %
                                     (anti_metabolite, fitness, (fitness_value - fitness)/fitness_value))

            new_fitness_values = {val: antimets for val, antimets in new_fitness_values.items()
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


class KnockoutMetaboliteFunction(object):
    def __init__(self, metabolites, ignore_transport=True, allow_accumulation=True):
        self.metabolites = metabolites
        self.ignore_transport = ignore_transport
        self.allow_accumulation = allow_accumulation

    def __call__(self, time_machine=None):
        for met in self.metabolites:
            knockout_metabolite(metabolite=met,
                                ignore_transport=self.ignore_transport,
                                allow_accumulation=self.allow_accumulation,
                                time_machine=time_machine)


def knockout_antimetabolite_substitution(model, reaction, ref_flux=0, ignore_metabolites=None, ignore_transport=True,
                                         allow_accumulation=True):
    """
    Generates a dictionary (species_id) -> function.
    The function blocks blocks the usage of a certain species by the model.

    Attributes
    ----------

    model: cobra.Model
    reaction: cobra.Reaction
    ref_flux: float
    ignore_metabolites: list
    ignore_transport: bool
    allow_accumulation: bool


    Returns
    -------
    dict

    """
    assert isinstance(reaction, Reaction)
    assert isinstance(ignore_metabolites, (list, set, tuple))

    if ref_flux != 0:
        substrates = [m for m, coeff in reaction.metabolites.items()
                      if coeff * ref_flux > 0 and m.id not in ignore_metabolites]
    else:
        if reaction.reversibility:
            substrates = [m for m in reaction.metabolites.keys() if m.id[:-2] not in ignore_metabolites]
        else:
            substrates = [m for m, coeff in reaction.metabolites.items()
                          if coeff > 0 and m.id[:-2] not in ignore_metabolites]

    species_ids = [m.id[:-2] for m in substrates]
    result = {}
    for species_id in species_ids:
        metabolites = search_metabolites(model, species_id)
        result[species_id] = KnockoutMetaboliteFunction(metabolites, ignore_transport, allow_accumulation)

    return result


def compete_metabolite(metabolite, reference_dist, fraction=0.5, ignore_transport=True,
                       allow_accumulation=True, time_machine=None):
    """
    Inhibits the usage of a metabolite based on a reference flux distributions.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A metabolite.
    reference_dist: dict or FluxDistributionResult
        The result of a FBA like simulation. Alternative can be dictionaries of reaction.id -> flux.
    fraction: float
        How much does it requires the reactions to go up.
    ignore_transport: bool
        Choose to ignore transport reactions.
    allow_accumulation: bool
        Allow to accumulate the metabolite (add a exchange reaction).
    time_machine: cameo.util.TimeMachine
        Action control.
    """
    reactions = metabolite.reactions

    if ignore_transport:
        reactions = [r for r in reactions if not len(set(m.compartment for m in r.metabolites)) > 1]

    if not isinstance(reference_dist, (dict, FluxDistributionResult)):
        raise ValueError("'reference_dist' must be a dict or FluxDistributionResult")

    exchanges = metabolite.model.exchanges

    for reaction in reactions:
        if reaction in exchanges:
            continue

        if reaction.metabolites[metabolite] > 0 and reference_dist[reaction.id] > 1e-6:
            reaction.change_bounds(lb=reference_dist[reaction.id] + (reference_dist[reaction.id] * fraction),
                                   time_machine=time_machine)

        if reaction.metabolites[metabolite] < 0 and reference_dist[reaction.id] < -1e-6:
            reaction.change_bounds(ub=reference_dist[reaction.id] + (reference_dist[reaction.id] * fraction),
                                   time_machine=time_machine)

    if allow_accumulation:
        species_id = metabolite.id[:-2]
        if "EX_%s_e" % species_id not in metabolite.model.reactions:
            return metabolite.model.add_exchange(metabolite, prefix="OF_", time_machine=time_machine)
        else:
            return metabolite.model.reactions.get_by_id("EX_%s_e" % species_id)


def inhibit_metabolite(metabolite, reference_dist, fraction=0.5, ignore_transport=True,
                       allow_accumulation=True, time_machine=None):
    """
    Inhibits the usage of a metabolite based on a reference flux distributions.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A metabolite.
    reference_dist: dict, FluxDistributionResult
        The result of a FBA like simulation. Alternative can be dictionaries of reaction.id -> flux.
    fraction: float, dict
        How much does it inhibits the reactions. A float applies the same amount of inhibition. A dictionary must
        contain an inhibition percentage to all reactions associated with the metabolite.
    ignore_transport: bool
        Choose to ignore transport reactions.
    allow_accumulation: bool
        Allow to accumulate the metabolite (add a exchange reaction).
    time_machine: cameo.util.TimeMachine
        Action control.
    """
    reactions = metabolite.reactions

    if ignore_transport:
        reactions = [r for r in reactions if not len(set(m.compartment for m in r.metabolites)) > 1]

    if isinstance(fraction, Number):
        fraction = {r.id: fraction for r in reactions}

    if not isinstance(fraction, dict):
        raise ValueError("'percentage' must be a dict or numeric")

    if not isinstance(reference_dist, (dict, FluxDistributionResult)):
        raise ValueError("'reference_dist' must be a dict or FluxDistributionResult")

    missing = []

    for reaction in reactions:
        if reaction.id not in fraction:
            missing.append(reaction)

    if len(missing):
        raise ValueError("reactions %s are not accounted in 'fraction'" % ", ".join("'%s'" % r.id for r in missing))

    exchanges = metabolite.model.exchanges

    for reaction in reactions:
        if reaction in exchanges:
            continue

        if reaction.metabolites[metabolite] < 0:
            if reference_dist[reaction.id] > 1e-6:
                reaction.change_bounds(ub=reference_dist[reaction.id] * fraction[reaction.id], time_machine=time_machine)
            else:
                reaction.change_bounds(ub=0, time_machine=time_machine)

        if reaction.metabolites[metabolite] > 0:
            if reference_dist[reaction.id] < -1e-6:
                reaction.change_bounds(lb=reference_dist[reaction.id] * fraction[reaction.id], time_machine=time_machine)
            else:
                reaction.change_bounds(lb=0, time_machine=time_machine)

    if allow_accumulation:
        species_id = metabolite.id[:-2]
        if "EX_%s_e" % species_id not in metabolite.model.reactions:
            return metabolite.model.add_exchange(metabolite, prefix="I_", demand=True, time_machine=time_machine)
        else:
            return metabolite.model.reactions.get_by_id("EX_%s_e" % species_id)


def knockout_metabolite(metabolite, ignore_transport=True, allow_accumulation=True, time_machine=None):
    """
    Inhibits the usage of a metabolite based on a reference flux distributions.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A metabolite.
    ignore_transport: bool
        Choose to ignore transport reactions.
    allow_accumulation: bool
        Allow to accumulate the metabolite (add a exchange reaction).
    time_machine: cameo.util.TimeMachine
        Action control.
    """
    reactions = metabolite.reactions

    if ignore_transport:
        reactions = [r for r in reactions if not len(set(m.compartment for m in r.metabolites)) > 1]

    exchanges = metabolite.model.exchanges

    for reaction in reactions:
        if reaction in exchanges:
            continue

        if reaction.reversibility:
            reaction.change_bounds(lb=0, ub=0, time_machine=time_machine)
        elif reaction.metabolites[metabolite] < 0:
            reaction.change_bounds(ub=0, time_machine=time_machine)

    if allow_accumulation:
        species_id = metabolite.id[:-2]
        if "EX_%s_e" % species_id not in metabolite.model.reactions:
            return metabolite.model.add_exchange(metabolite, prefix="KO_", demand=True, time_machine=time_machine)
        else:
            return metabolite.model.reactions.get_by_id("EX_%s_e" % species_id)