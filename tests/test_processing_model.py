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
import unittest

from cameo import load_model, pfba
from cameo.util import TimeMachine

from marsi.processing.models import degree, in_degree, out_degree
from marsi.processing.models import inhibit_metabolite, compete_metabolite, knockout_metabolite
from marsi.processing.models import apply_antimetabolite, knockout_antimetabolite_substitution


TEST_DIR = os.path.dirname(__file__)


class ModelProcessingTestCase(unittest.TestCase):
    iaf_model = None

    @classmethod
    def setUpClass(cls):
        cls.iaf_model = load_model(os.path.join(TEST_DIR, 'fixtures', 'iAF1260.json'))
        cls.iaf_reference = pfba(cls.iaf_model)
        cls.ijo_model = load_model(os.path.join(TEST_DIR, 'fixtures', 'iJO1366.json'))
        cls.ijo_reference = pfba(cls.ijo_model)

    def metabolite_degree_test(self):
        _2ddecg3p_c = self.iaf_model.metabolites.get_by_id('2ddecg3p_c')
        self.assertEqual(in_degree(_2ddecg3p_c), 1)
        self.assertEqual(out_degree(_2ddecg3p_c), 1)
        self.assertEqual(degree(_2ddecg3p_c), 2)

        g3p_c = self.iaf_model.metabolites.g3p_c
        self.assertEqual(in_degree(g3p_c), 10)
        self.assertEqual(out_degree(g3p_c), 13)
        self.assertEqual(degree(g3p_c), 23)

    def inhibit_metabolite_test(self):
        succ_c = self.iaf_model.metabolites.succ_c

    def knockout_metabolite_with_exchange_test(self):
        succ_c = self.iaf_model.metabolites.succ_c
        succ_e = self.iaf_model.metabolites.succ_e
        succ_p = self.iaf_model.metabolites.succ_p

        succ_c_transport_in = ["SUCCt2_2pp", "SUCCt2_3pp"]
        succ_c_transport_out = ["CITt7pp",  "SUCCt3pp"]
        succ_c_transport_both = ["SUCFUMtpp", "TARTRt7pp"]

        succ_c_transport_reactions = [self.iaf_model.reactions.get_by_id(rid) for rid in succ_c_transport_in]
        succ_c_transport_reactions += [self.iaf_model.reactions.get_by_id(rid) for rid in succ_c_transport_out]
        succ_c_transport_reactions += [self.iaf_model.reactions.get_by_id(rid) for rid in succ_c_transport_both]
        succ_c_transport_bounds = {r.id: (r.lower_bound, r.upper_bound) for r in succ_c_transport_reactions}

        with TimeMachine() as tm:
            knockout_metabolite(succ_c, ignore_transport=True, allow_accumulation=True, time_machine=tm)
            self.assertNotIn("KO_succ_c", self.iaf_model.reactions)
            for transport_r in succ_c_transport_reactions:
                self.assertEqual(transport_r.lower_bound, succ_c_transport_bounds[transport_r.id][0])
                self.assertEqual(transport_r.upper_bound, succ_c_transport_bounds[transport_r.id][1])

        with TimeMachine() as tm:
            knockout_metabolite(succ_c, ignore_transport=True, allow_accumulation=False, time_machine=tm)
            self.assertNotIn("KO_succ_c", self.iaf_model.reactions)
            for transport_r in succ_c_transport_reactions:
                self.assertEqual(transport_r.lower_bound, succ_c_transport_bounds[transport_r.id][0])
                self.assertEqual(transport_r.upper_bound, succ_c_transport_bounds[transport_r.id][1])

        with TimeMachine() as tm:
            knockout_metabolite(succ_c, ignore_transport=False, allow_accumulation=True, time_machine=tm)
            self.assertNotIn("KO_succ_c", self.iaf_model.reactions)
            for transport_rid in succ_c_transport_out:
                transport_r = self.iaf_model.reactions.get_by_id(transport_rid)
                self.assertEqual(transport_r.lower_bound, succ_c_transport_bounds[transport_r.id][0])
                self.assertEqual(transport_r.upper_bound, 0)

            for transport_rid in succ_c_transport_both:
                transport_r = self.iaf_model.reactions.get_by_id(transport_rid)
                self.assertEqual(transport_r.lower_bound, 0)
                self.assertEqual(transport_r.upper_bound, 0)

        with TimeMachine() as tm:
            knockout_metabolite(succ_e, ignore_transport=True, allow_accumulation=True, time_machine=tm)
            self.assertNotIn("KO_succ_p", self.iaf_model.reactions)

        with TimeMachine() as tm:
            knockout_metabolite(succ_e, ignore_transport=False, allow_accumulation=True, time_machine=tm)
            self.assertNotIn("KO_succ_p", self.iaf_model.reactions)
            for transport_rid in succ_c_transport_in:
                transport_r = self.iaf_model.reactions.get_by_id(transport_rid)
                self.assertEqual(transport_r.lower_bound, succ_c_transport_bounds[transport_r.id][0])
                self.assertEqual(transport_r.upper_bound, 0)

            for transport_rid in succ_c_transport_both:
                transport_r = self.iaf_model.reactions.get_by_id(transport_rid)
                self.assertEqual(transport_r.lower_bound, 0)
                self.assertEqual(transport_r.upper_bound, 0)

        with TimeMachine() as tm:
            knockout_metabolite(succ_p, ignore_transport=True, allow_accumulation=True, time_machine=tm)
            self.assertNotIn("KO_succ_p", self.iaf_model.reactions)

    def knockout_metabolite_without_exchange_test(self):
        raise NotImplementedError

    def compete_metabolite_test(self):
        ala__L_c = self.iaf_model.metabolites.ala__L_c

    def knockout_antimetabolite_substitution_test(self):
        targets = []

