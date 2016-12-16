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

import pybel
from marsi.io.mongodb import Metabolite
from marsi.processing.chemistry.openbabel import fingerprint_to_bits


class MetaboliteTestCase(unittest.TestCase):
    def setUp(self):
        from mongoengine import connect
        connect("marsi-db")

    def test_fingerprint(self):
        met1 = Metabolite.objects[0]
        fp = met1.fingerprint('maccs')
        ob_fp = met1.molecule(library='openbabel').calcfp('maccs')
        self.assertTrue((fp == ob_fp.fp).all())

        fp = met1._fingerprints['maccs']
        uint_vector = pybel.ob.vectorUnsignedInt(len(fp))
        for i in range(len(fp)):
            uint_vector[i] = fp[i]
        regenerated_fp = pybel.Fingerprint(uint_vector)

        self.assertEqual(regenerated_fp | ob_fp, 1.0)

        bits = met1.fingerprint('maccs', hash=True, bits=165)
        ob_bits = fingerprint_to_bits(ob_fp, bits=165)
        for i in range(165):
            self.assertEqual(bits[i], ob_bits[i], msg="bit %i" % i)
