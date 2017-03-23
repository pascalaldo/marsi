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
import pytest

from marsi.cobra.strain_design.evolutionary import OptMet, OptMetResult

CURRENT_DIRECTORY = os.path.dirname(__file__)
FIXTURES = os.path.join(CURRENT_DIRECTORY, 'fixtures')


def test_succinate(model):
    optimization = OptMet(model=model, plot=False)

    optimization_kwargs = dict(max_evaluations=1000, max_knockouts=6, target="succ_e",
                               substrate="EX_glc__D_e", biomass=model.biomass)

    assert optimization.manipulation_type == "metabolites"

    # result = optimization.run(**optimization_kwargs)
    #
    # assert isinstance(result, OptMetResult)
    # assert len(result) > 0
