# Copyright 2017 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
from numbers import Number

from cameo.flux_analysis.simulation import FluxDistributionResult


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
