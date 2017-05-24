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
from cameo import fba
from cameo.core.strain_design import StrainDesign
from cameo.core.target import ReactionKnockoutTarget, ReactionModulationTarget
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
from cameo.util import TimeMachine, flatten

from marsi.cobra.strain_design.post_processing import find_anti_metabolite_knockouts, find_anti_metabolite_modulation, \
    convert_target, replace_design
from marsi.cobra.strain_design.target import MetaboliteKnockoutTarget, AntiMetaboliteManipulationTarget


def test_find_anti_metabolite_knockouts(model, essential_metabolites):
    gapd = model.reactions.GAPD

    anti_metabolites = find_anti_metabolite_knockouts(gapd, 0.5, essential_metabolites)
    assert all(isinstance(anti_met, MetaboliteKnockoutTarget) for anti_met in anti_metabolites.values())


def test_find_anti_metabolite_modulation(model, essential_metabolites):
    gapd = model.reactions.GAPD

    anti_metabolites = find_anti_metabolite_modulation(gapd, 0.5, essential_metabolites)
    assert all(isinstance(anti_met, AntiMetaboliteManipulationTarget) for anti_met in anti_metabolites.values())


def test_convert_target(model, essential_metabolites):
    ko_target = ReactionKnockoutTarget("ACALD")
    mod_target = ReactionModulationTarget("ACALD", 0.1, 0.4)

    converted_ko_targets = convert_target(model, ko_target, essential_metabolites)
    assert all(isinstance(anti_met, MetaboliteKnockoutTarget) for anti_met in converted_ko_targets.values())

    converted_mod_targets = convert_target(model, mod_target, essential_metabolites)
    assert all(isinstance(anti_met, AntiMetaboliteManipulationTarget) for anti_met in converted_mod_targets.values())


def test_convert_design(model, essential_metabolites):
    # Target: EX_lac__D_e
    # Medium: glucose
    # Aerobicity: anaerobic
    # fva_min: 16.8641
    # Design: ΔACALD ΔALCD2x ΔPTA2 ΔPTAr

    targets = [
        ReactionKnockoutTarget("ACALD"),
        ReactionKnockoutTarget("ALCD2x"),
        ReactionKnockoutTarget("PTA2"),
        ReactionKnockoutTarget("PTAr")
    ]

    strain_design = StrainDesign(targets)
    objective_function = biomass_product_coupled_yield(model.biomass, "EX_lac__D_e", "EX_glc__D_e")
    model.reactions.EX_o2_e.lower_bound = 0

    with TimeMachine() as tm:
        strain_design.apply(model, time_machine=tm)
        solution = fba(model, objective=model.biomass)
        fitness = objective_function(model, solution, targets)

    replacement = replace_design(model, strain_design, fitness, objective_function,
                                 fba, essential_metabolites=essential_metabolites)

    print(replacement)
    model.reactions.EX_o2_e.lower_bound = 1000

    assert "actp" in [target.id for target in flatten(replacement.metabolite_targets.values)]
