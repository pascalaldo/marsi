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

import pytest
from cameo import pfba, fba
from cameo.util import TimeMachine

from marsi.cobra.flux_analysis.analysis import sensitivity_analysis
from marsi.utils import search_metabolites
from marsi.cobra.flux_analysis.manipulation import knockout_metabolite, compete_metabolite, inhibit_metabolite


@pytest.fixture(params=[True, False])
def allow_accumulation(request):
    return request.param


@pytest.fixture(params=[True, False])
def ignore_transport(request):
    return request.param


def test_search_metabolites(model):
    results = search_metabolites(model, "glc__D", ignore_external=True)
    assert any(met.id[-2:] != "_e" for met in results)

    results = search_metabolites(model, "glc__D", ignore_external=False)
    assert any(met.id[-2:] == "_e" for met in results)


def test_inhibit_metabolite(model, allow_accumulation, benchmark):
    succ_c = model.metabolites.succ_c

    reference = pfba(model, objective=model.biomass)
    time_machine = TimeMachine()

    def _inhibit_metabolite(tm):
        inhibit_metabolite(model, succ_c, reference, allow_accumulation=allow_accumulation, time_machine=tm)
        tm.reset()

    benchmark(_inhibit_metabolite, time_machine)
    with time_machine as tm:
        exchange = inhibit_metabolite(model, succ_c, reference, allow_accumulation=allow_accumulation, time_machine=tm)

        result = pfba(model, objective=model.biomass)

    reference_consumption_turnover = 0
    result_consumption_turnover = 0
    for reaction in succ_c.reactions:
        coefficient = reaction.metabolites[succ_c]
        ref_t = coefficient * reference[reaction]
        res_t = coefficient * result[reaction]

        if ref_t < 0:
            reference_consumption_turnover += ref_t

        if res_t < 0:
            result_consumption_turnover += res_t

    print(reference_consumption_turnover, result_consumption_turnover)
    if allow_accumulation:
        print(exchange.id, result[exchange])
    assert abs(reference_consumption_turnover) > abs(result_consumption_turnover)


def test_sensitivity_analysis(model, benchmark):
    succ_c = model.metabolites.succ_c

    reference = pfba(model, objective=model.biomass)

    result = benchmark.pedantic(sensitivity_analysis,
                                args=(model, succ_c, model.biomass),
                                kwargs=dict(is_essential=False, reference=reference),
                                iterations=1,
                                rounds=1)

    print(result)


def test_knockout_metabolite_knockout_exchangeable(model, allow_accumulation, ignore_transport, benchmark):
    succ_c = model.metabolites.succ_c

    succ_c_transport_in = ["SUCCt2_2pp", "SUCCt2_3pp"]
    succ_c_transport_out = ["CITt7pp",  "SUCCt3pp"]
    succ_c_transport_both = ["SUCFUMtpp", "TARTRt7pp"]

    succ_c_transport_reactions = [model.reactions.get_by_id(rid) for rid in succ_c_transport_in]
    succ_c_transport_reactions += [model.reactions.get_by_id(rid) for rid in succ_c_transport_out]
    succ_c_transport_reactions += [model.reactions.get_by_id(rid) for rid in succ_c_transport_both]
    succ_c_transport_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in succ_c_transport_reactions}

    with TimeMachine() as tm:

        benchmark(knockout_metabolite, model, succ_c,
                  ignore_transport=ignore_transport,
                  allow_accumulation=allow_accumulation,
                  time_machine=tm)

        for transport_r in succ_c_transport_reactions:
            if ignore_transport:
                assert transport_r.lower_bound == succ_c_transport_bounds[transport_r.id][0]
                assert transport_r.upper_bound == succ_c_transport_bounds[transport_r.id][1]
            else:
                if transport_r in succ_c_transport_in:
                    assert transport_r.lower_bound == 0
                    assert transport_r.upper_bound == succ_c_transport_bounds[transport_r.id][1]
                elif transport_r in succ_c_transport_out:
                    assert transport_r.lower_bound == succ_c_transport_bounds[transport_r.id][1]
                    assert transport_r.upper_bound == 0
                elif transport_r in succ_c_transport_both:
                    assert transport_r.lower_bound == 0
                    assert transport_r.upper_bound == 0


def test_knockout_metabolite_knockout_non_exchangeable(model, allow_accumulation, ignore_transport, benchmark):
    xu5p__L = model.metabolites.xu5p__L_c

    def _knockout_metabolite():
        with TimeMachine() as tm:
            knockout_metabolite(model, xu5p__L, ignore_transport=ignore_transport,
                                allow_accumulation=allow_accumulation, time_machine=tm)

    benchmark(_knockout_metabolite)

    with TimeMachine() as tm:
        knockout_metabolite(model, xu5p__L, ignore_transport=ignore_transport,
                            allow_accumulation=allow_accumulation, time_machine=tm)

        if allow_accumulation:
            assert "KO_xu5p__L_c" in model.reactions


def test_compete_metabolite_test(model, amino_acid, benchmark):
    aa = model.metabolites.get_by_id(amino_acid)
    reference = pfba(model, objective=model.biomass)

    time_machine = TimeMachine()

    def _compete_metabolite(time_machine):
        compete_metabolite(model, aa, reference, time_machine=time_machine)
        time_machine.reset()

    try:
        benchmark(_compete_metabolite, time_machine)
    except ValueError:
        assert fba(model, objective=aa).objective_value <= 1e-6
        return

    with time_machine as tm:
        exchange = compete_metabolite(model, aa, fraction=0.1, reference_dist=reference, time_machine=tm)

        solution = pfba(model, objective=model.biomass)

    assert solution[exchange] > 0
