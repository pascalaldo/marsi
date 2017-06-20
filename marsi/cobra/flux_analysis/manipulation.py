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

import logging

from cameo.util import TimeMachine
from cobra import DictList
from functools import partial

from cameo.flux_analysis.simulation import FluxDistributionResult

__all__ = ["compete_metabolite", "inhibit_metabolite", "knockout_metabolite", "apply_anti_metabolite"]

logger = logging.getLogger(__name__)


def compete_metabolite(model, metabolite, reference_dist, fraction=0.5, allow_accumulation=True,
                       constant=1e4, time_machine=None):
    """
    Increases the usage of a metabolite based on a reference flux distributions.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        A constraint-based model.
    metabolite : cobra.Metabolite
        A metabolite.
    reference_dist : dict or FluxDistributionResult
        The result of a FBA like simulation. Alternative can be dictionaries of reaction.id -> flux.
    fraction : float
        How much does it requires the reactions to go up.
    allow_accumulation : bool
        Allow to accumulate the metabolite (add a exchange reaction).
    constant : float
        A large number (like 10000).
    time_machine : cameo.util.TimeMachine
        Action control.

    Returns
    -------
    cameo.core.Reaction
        If allow accumulation returns the exchange reaction associated with the metabolite.
    """

    reactions = [r for r in metabolite.reactions if len(set(m.compartment for m in r.metabolites)) == 1]

    if isinstance(reference_dist, FluxDistributionResult):
        reference_dist = reference_dist.fluxes

    if not isinstance(reference_dist, dict):
        raise ValueError("'reference_dist' must be a dict or FluxDistributionResult")

    exchanges = DictList(model.exchanges)

    if allow_accumulation:
        species_id = metabolite.id[:-2]
        if "EX_%s_e" % species_id not in exchanges:
            exchange = model.add_exchange(metabolite, prefix="COMPETE_", time_machine=time_machine)
        else:
            exchange = model.reactions.get_by_id("EX_%s_e" % species_id)
    else:
        exchange = None

    aux_variables = {}
    ind_variables = {}
    turnover = sum(abs(r.metabolites[metabolite] * reference_dist.get(r.id, 0)) for r in metabolite.reactions)
    for reaction in reactions:
        coefficient = reaction.metabolites[metabolite]

        # Optimization to reduce y variables and problem complexity:
        # Irreversible reactions that only produce the metabolite can be ignored because they will not contribute
        # to the consumption turnover. Reactions that only consume the metabolite can be added directly into the
        # sum constraint. This allows for a less complex problem with less variables.

        if not reaction.reversibility:
            if coefficient < 0:  # skip reactions that can only produce the metabolite
                continue
            else:  # keep the v*coefficient value for reactions that can only consume the metabolite
                aux_variables[reaction.id] = reaction.flux_expression * coefficient
                continue

        to_add = []

        ind_var_id = "y_%s" % reaction.id
        aux_var_id = "u_%s" % reaction.id
        try:
            ind_var = model.solver.variables[ind_var_id]
            aux_var = model.solver.variables[aux_var_id]
        except KeyError:
            ind_var = model.solver.interface.Variable(ind_var_id, type='binary')
            aux_var = model.solver.interface.Variable(aux_var_id, lb=0)
            to_add += [ind_var, aux_var]

        aux_variables[reaction.id] = aux_var
        ind_variables[reaction.id] = ind_var

        upper_indicator_constraint_name = "ind_%s_u" % reaction.id
        lower_indicator_constraint_name = "ind_%s_l" % reaction.id

        auxiliary_constraint_a_name = "aux_%s_a" % reaction.id
        auxiliary_constraint_b_name = "aux_%s_b" % reaction.id
        auxiliary_constraint_c_name = "aux_%s_c" % reaction.id
        auxiliary_constraint_d_name = "aux_%s_d" % reaction.id

        try:
            model.solver.constraints[upper_indicator_constraint_name]
        except KeyError:

            # constraint y to be 0 if Sv >= 0 (production)

            #   -M             0                M
            # Sv <-------------|---------------->
            #         y=0      |      y=1

            # Sv - My <= 0
            # if y = 1 then Sv <= M
            # if y = 0 then Sv <= 0
            upper_indicator_expression = coefficient * reaction.flux_expression - ind_var * constant
            ind_constraint_u = model.solver.interface.Constraint(upper_indicator_expression,
                                                                 name=upper_indicator_constraint_name,
                                                                 ub=0)

            # Sv + M(1-y) >= 0
            # if y = 1 then Sv >= 0
            # if y = 0 then Sv >= -M
            lower_indicator_expression = coefficient * reaction.flux_expression + constant - ind_var * constant
            ind_constraint_l = model.solver.interface.Constraint(lower_indicator_expression,
                                                                 name=lower_indicator_constraint_name,
                                                                 lb=0)

            # a) -My + u <= 0
            # b) My + u >= 0

            # if y = 0, u = 0
            # if y = 1, -M <= u <= M
            aux_indicator_expression_a = -constant * ind_var + aux_var
            aux_constraint_a = model.solver.interface.Constraint(aux_indicator_expression_a,
                                                                 name=auxiliary_constraint_a_name,
                                                                 ub=0)

            aux_indicator_expression_b = constant * ind_var + aux_var
            aux_constraint_b = model.solver.interface.Constraint(aux_indicator_expression_b,
                                                                 name=auxiliary_constraint_b_name,
                                                                 lb=0)
            #
            # # c) -M(1-y) + u - viSi <= 0
            # # d) M(1-y) + u - viSi >= 0
            #
            # # if y = 1 then 0 <= u - viSi <= 0
            # # if y = 0 then -M <= u - viSi <= M
            aux_indicator_expression_c = -constant * (1 - ind_var) + aux_var - reaction.flux_expression * coefficient
            aux_constraint_c = model.solver.interface.Constraint(aux_indicator_expression_c,
                                                                 name=auxiliary_constraint_c_name,
                                                                 ub=0)

            aux_indicator_expression_d = constant * (1 - ind_var) + aux_var - reaction.flux_expression * coefficient
            aux_constraint_d = model.solver.interface.Constraint(aux_indicator_expression_d,
                                                                 name=auxiliary_constraint_d_name,
                                                                 lb=0)

            to_add += [ind_constraint_l, ind_constraint_u, aux_constraint_a,
                       aux_constraint_b, aux_constraint_c, aux_constraint_d]

            if time_machine:
                time_machine(do=partial(model.solver.add, to_add), undo=partial(model.solver.remove, to_add))
            else:
                model.solver.add(to_add)

    model.solver.update()
    min_production_turnover = (1 + fraction) * (turnover / 2)
    # sum(u) >= (1 + fraction) * uWT
    decrease_turnover_constraint = model.solver.interface.Constraint(sum(aux_variables.values()),
                                                                     name="take_less_%s" % metabolite.id,
                                                                     lb=min_production_turnover)

    if time_machine:
        time_machine(do=partial(model.solver.add, decrease_turnover_constraint),
                     undo=partial(model.solver.remove, decrease_turnover_constraint))
    else:
        model.solver.add(decrease_turnover_constraint)

    return exchange


def inhibit_metabolite(model, metabolite, reference_dist, fraction=0.5, allow_accumulation=True,
                       constant=1e4, time_machine=None):
    """
    Inhibits the usage of a metabolite based on a reference flux distributions.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        A constraint-based model.
    metabolite : cobra.Metabolite
        A metabolite.
    reference_dist : dict, FluxDistributionResult
        The result of a FBA like simulation. Alternative can be dictionaries of reaction.id -> flux.
    fraction : float
        How much does it inhibits the reactions. A float applies the same amount of inhibition. A dictionary must
        contain an inhibition percentage to all reactions associated with the metabolite.
    allow_accumulation : bool
        Allow to accumulate the metabolite (add a exchange reaction).
    constant : float
        A large number (like 10000).
    time_machine : cameo.util.TimeMachine
        Action control.

    Returns
    -------
    cameo.core.Reaction, None
        If allow accumulation returns the exchange reaction associated with the metabolite.

    """
    reactions = [r for r in metabolite.reactions if len(set(m.compartment for m in r.metabolites)) == 1]

    if isinstance(reference_dist, FluxDistributionResult):
        reference_dist = reference_dist.fluxes

    if not isinstance(reference_dist, dict):
        raise ValueError("'reference_dist' must be a dict or FluxDistributionResult")

    exchanges = DictList(model.exchanges)

    if allow_accumulation:
        species_id = metabolite.id[:-2]
        if "EX_%s_e" % species_id not in exchanges:
            exchange = model.add_exchange(metabolite, prefix="INHIBIT_", time_machine=time_machine)
        else:
            exchange = model.reactions.get_by_id("EX_%s_e" % species_id)
    else:
        exchange = None

    aux_variables = {}
    ind_variables = {}
    turnover = sum(abs(r.metabolites[metabolite] * reference_dist.get(r.id, 0)) for r in metabolite.reactions)
    for reaction in reactions:
        coefficient = reaction.metabolites[metabolite]

        # Optimization to reduce y variables and problem complexity:
        # Irreversible reactions that only produce the metabolite can be ignored because they will not contribute
        # to the consumption turnover. Reactions that only consume the metabolite can be added directly into the
        # sum constraint. This allows for a less complex problem with less variables.

        if not reaction.reversibility:
            if coefficient > 0:  # skip reactions that can only produce the metabolite
                continue
            else:  # keep the v*coefficient value for reactions that can only consume the metabolite
                aux_variables[reaction.id] = - reaction.flux_expression * coefficient
                continue

        to_add = []

        ind_var_id = "y_%s" % reaction.id
        aux_var_id = "u_%s" % reaction.id
        try:
            ind_var = model.solver.variables[ind_var_id]
            aux_var = model.solver.variables[aux_var_id]
        except KeyError:
            ind_var = model.solver.interface.Variable(ind_var_id, type='binary')
            aux_var = model.solver.interface.Variable(aux_var_id, lb=0)
            to_add += [ind_var, aux_var]

        aux_variables[reaction.id] = aux_var
        ind_variables[reaction.id] = ind_var

        upper_indicator_constraint_name = "ind_%s_u" % reaction.id
        lower_indicator_constraint_name = "ind_%s_l" % reaction.id

        auxiliary_constraint_a_name = "aux_%s_a" % reaction.id
        auxiliary_constraint_b_name = "aux_%s_b" % reaction.id
        auxiliary_constraint_c_name = "aux_%s_c" % reaction.id
        auxiliary_constraint_d_name = "aux_%s_d" % reaction.id

        try:
            model.solver.constraints[upper_indicator_constraint_name]
        except KeyError:

            # constraint y to be 0 if Sv >= 0 (production)

            #   -M             0                M
            # Sv <-------------|---------------->
            #         y=0      |      y=1

            # -Sv - My <= 0
            # if y = 1 then Sv <= M
            # if y = 0 then Sv >= 0
            upper_indicator_expression = - coefficient * reaction.flux_expression - ind_var * constant
            ind_constraint_u = model.solver.interface.Constraint(upper_indicator_expression,
                                                                 name=upper_indicator_constraint_name,
                                                                 ub=0)

            # -Sv + M(1-y) >= 0
            # if y = 1 then Sv <= 0
            # if y = 0 then Sv <= M
            lower_indicator_expression = - coefficient * reaction.flux_expression + constant - ind_var * constant
            ind_constraint_l = model.solver.interface.Constraint(lower_indicator_expression,
                                                                 name=lower_indicator_constraint_name,
                                                                 lb=0)

            # a) -My + u <= 0
            # b) My + u >= 0

            # if y = 0, u = 0
            # if y = 1, -M <= u <= M
            aux_indicator_expression_a = -constant * ind_var + aux_var
            aux_constraint_a = model.solver.interface.Constraint(aux_indicator_expression_a,
                                                                 name=auxiliary_constraint_a_name,
                                                                 ub=0)

            aux_indicator_expression_b = constant * ind_var + aux_var
            aux_constraint_b = model.solver.interface.Constraint(aux_indicator_expression_b,
                                                                 name=auxiliary_constraint_b_name,
                                                                 lb=0)
            #
            # # c) -M(1-y) + u + viSi <= 0
            # # d) M(1-y) + u + viSi >= 0
            #
            # # if y = 1 then 0 <= u + viSi <= 0
            # # if y = 0 then -M <= u + viSi <= M
            aux_indicator_expression_c = -constant * (1 - ind_var) + aux_var + reaction.flux_expression * coefficient
            aux_constraint_c = model.solver.interface.Constraint(aux_indicator_expression_c,
                                                                 name=auxiliary_constraint_c_name,
                                                                 ub=0)

            aux_indicator_expression_d = constant * (1 - ind_var) + aux_var + reaction.flux_expression * coefficient
            aux_constraint_d = model.solver.interface.Constraint(aux_indicator_expression_d,
                                                                 name=auxiliary_constraint_d_name,
                                                                 lb=0)

            to_add += [ind_constraint_l, ind_constraint_u, aux_constraint_a,
                       aux_constraint_b, aux_constraint_c, aux_constraint_d]

            if time_machine:
                time_machine(do=partial(model.solver.add, to_add), undo=partial(model.solver.remove, to_add))
            else:
                model.solver.add(to_add)

    model.solver.update()
    max_production_turnover = (1 - fraction) * (turnover / 2)
    # sum(u) <= (1-fraction) * uWT
    decrease_turnover_constraint = model.solver.interface.Constraint(sum(aux_variables.values()),
                                                                     name="take_less_%s" % metabolite.id,
                                                                     ub=max_production_turnover)

    if time_machine:
        time_machine(do=partial(model.solver.add, decrease_turnover_constraint),
                     undo=partial(model.solver.remove, decrease_turnover_constraint))
    else:
        model.solver.add(decrease_turnover_constraint)

    return exchange


def knockout_metabolite(model, metabolite, ignore_transport=True, allow_accumulation=True, time_machine=None):
    """
    Inhibits the usage of a metabolite based on a reference flux distributions.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        A constraint-based model.
    metabolite: cobra.Metabolite
        A metabolite.
    ignore_transport : bool
        Choose to ignore transport reactions.
    allow_accumulation : bool
        Allow to accumulate the metabolite (add a exchange reaction).
    time_machine : cameo.util.TimeMachine
        Action control.

    Returns
    -------
    cameo.core.Reaction, None
        If allow accumulation returns the exchange reaction associated with the metabolite.
    """
    reactions = metabolite.reactions

    if ignore_transport:
        reactions = [r for r in reactions if not len(set(m.compartment for m in r.metabolites)) > 1]

    exchanges = model.exchanges

    for reaction in reactions:
        if reaction in exchanges:
            continue

        if reaction.reversibility:
            reaction.change_bounds(lb=0, ub=0, time_machine=time_machine)
        elif reaction.metabolites[metabolite] < 0:
            reaction.change_bounds(ub=0, time_machine=time_machine)

    if allow_accumulation:
        species_id = metabolite.id[:-2]
        if "EX_%s_e" % species_id not in model.reactions:
            return model.add_exchange(metabolite, prefix="KO_", demand=True, time_machine=time_machine)
        else:
            return model.reactions.get_by_id("EX_%s_e" % species_id)


def apply_anti_metabolite(model, metabolites, essential_metabolites, reference, inhibition_fraction=.0,
                          competition_fraction=.0, allow_accumulation=True, time_machine=None):
    """
    Apply a metabolite in the context of a model without knowing if it is activating or inhibiting.

    Parameters
    ----------
    model : cameo.core.SolverBasedModel
        A constraint-based model.
    metabolites : list
        Metabolites of the same species.
    essential_metabolites : list
        A list of essential metabolites.
    reference : dict, cameo.core.FluxDistributionResult
        A flux distribution.
    inhibition_fraction : float
        How much a metabolite inhibits.
    competition_fraction : float
        How much a metabolite competes.
    allow_accumulation : bool
        Allow accumulation of the metabolite.
    time_machine : TimeMachine
        A context manager.

    Returns
    -------
    set
        Exchange reactions added for accumulation.

    """
    exchanges = set()

    if any(met in essential_metabolites for met in metabolites):
        for metabolite in metabolites:
            exchanges.add(compete_metabolite(model,
                                             metabolite,
                                             reference,
                                             allow_accumulation=allow_accumulation,
                                             fraction=competition_fraction,
                                             time_machine=time_machine))
    else:
        for metabolite in metabolites:
            exchanges.add(inhibit_metabolite(model,
                                             metabolite,
                                             reference,
                                             allow_accumulation=allow_accumulation,
                                             fraction=inhibition_fraction,
                                             time_machine=time_machine))

    return exchanges
