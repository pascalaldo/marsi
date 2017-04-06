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

import os
from cameo import fba
from cameo.core.strain_design import StrainDesign
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield


from marsi.cobra.strain_design.evolutionary import OptMet, OptMetResult, process_metabolite_knockout_solution

CURRENT_DIRECTORY = os.path.dirname(__file__)
FIXTURES = os.path.join(CURRENT_DIRECTORY, 'fixtures')


def test_design_processing_function(model):
    orignal_oxigen_uptake = model.reactions.EX_o2_e.lower_bound
    target = "EX_succ_e"
    substrate = "EX_glc__D_e"
    objective_function = biomass_product_coupled_yield(model.biomass, target, substrate)
    solution = ["mal__D"]
    try:
        model.reactions.EX_o2_e.lower_bound = 0
        result = process_metabolite_knockout_solution(model, solution, fba, {}, model.biomass, target,
                                                      substrate, objective_function)
    finally:
        model.reactions.EX_o2_e.lower_bound = orignal_oxigen_uptake

    design, size, fva_min, fva_max, target_flux, biomass_flux, _yield, fitness = result

    assert isinstance(design, StrainDesign)
    assert size == len(solution)
    assert size == 1
    assert fva_min > 0
    assert fva_max >= fva_min
    assert target_flux > 0
    assert biomass_flux > 0
    assert _yield > 0
    assert fitness > 0


def test_succinate(model):
    optimization = OptMet(model=model, plot=False)

    optimization_kwargs = dict(max_evaluations=1000, max_knockouts=6, target="succ_e",
                               substrate="EX_glc__D_e", biomass=model.biomass)

    assert optimization.manipulation_type == "metabolites"

    # result = optimization.run(**optimization_kwargs)
    #
    # assert isinstance(result, OptMetResult)
    # assert len(result) > 0
