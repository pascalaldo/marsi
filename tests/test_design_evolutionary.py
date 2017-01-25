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

import unittest

from cameo import models

from marsi.cobra.strain_design.evolutionary import OptGeneMet


class OptGeneMetTestCase(unittest.TestCase):
    model = None

    @classmethod
    def setUpClass(cls):
        cls.model = models.bigg.iJO1366

    def test_succinate(self):
        optimization = OptGeneMet(model=self.model, plot=False)
        result = optimization.run(max_evaluations=1000,
                                  max_knockouts=6,
                                  target="succ_e",
                                  substrate="EX_glc__D_e",
                                  biomass="BIOMASS_Ec_iJO1366_core_53p95M")

        self.assertTrue(any('succ' in metabolites for metabolites in result.data_frame.metabolites))
