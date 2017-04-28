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

import os
import pytest

from cameo import load_model

TEST_DIR = os.path.dirname(__file__)


BIOMASS_IDS = {
    "iJO1366": "BIOMASS_Ec_iJO1366_core_53p95M",
    "iAF1260": "BIOMASS_Ec_iAF1260_core_59p81M"
}


@pytest.fixture(params=["iJO1366", "iAF1260"], scope="global")
def model(request):
    """
    Genome-scale metabolic model. Loaded using cameo.load_model.

    Returns
    -------
    cameo.SolverBasedModel
    """
    model = load_model(os.path.join(TEST_DIR, 'fixtures', '%s.json' % request.param))
    setattr(model, 'biomass', BIOMASS_IDS[request.param])

    return model


@pytest.fixture(params=["ala__L_c", "ser__L_c", "trp__L_c", "glu__L_c"], scope='global')
def amino_acid(request):
    return request.param


@pytest.fixture(scope='global')
def essential_metabolites(model):
    return model.essential_metabolites()
