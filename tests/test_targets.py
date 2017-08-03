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

import pytest
from cameo import pfba

from marsi.cobra.strain_design.target import AntiMetaboliteManipulationTarget, MetaboliteKnockoutTarget


@pytest.fixture(params=["glc__D", "try__L", "atp"])
def species(request):
    return request.param


def test_anti_metabolite_manipulation_target(model, species):
    reference = pfba(model, objective=model.biomass)
    target = AntiMetaboliteManipulationTarget(species)
    compartments = model.compartments
    expected_ids = [species + "_" + compartment for compartment in compartments]
    metabolites = target.metabolites(model)

    assert all(m.id in expected_ids for m in metabolites)
    assert repr(target) == "<AntiMetaboliteManipulation %s (%.3f)>" % (target.id, target.fraction)
    assert str(target) == b'\xe2\x98\xa3'.decode('utf-8') + "(%.3f)-%s" % (target.fraction, target.id)

    with model:
        target.apply(model, reference)


def test_metabolite_knockout_target(model, species):
    target = MetaboliteKnockoutTarget(species)
    compartments = model.compartments
    expected_ids = [species + "_" + compartment for compartment in compartments]
    metabolites = target.metabolites(model)

    assert target.fraction == 0
    assert all(m.id in expected_ids for m in metabolites)
    assert repr(target) == "<MetaboliteKnockout %s>" % target.id
    assert str(target) == b'\xe2\x98\xa3'.decode('utf-8') + "-%s" % target.id

    with model:
        target.apply(model)
