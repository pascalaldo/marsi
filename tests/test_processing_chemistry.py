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

import numpy as np

from marsi.chemistry.common import SOLUBILITY, tanimoto_coefficient, tanimoto_distance
from marsi.chemistry import openbabel, rdkit

CARBON_ATOMIC_NUMBER = 6
HYDROGEN_ATOMIC_NUMBER = 1


class SolubilityTestCase(unittest.TestCase):
    def setUp(self):
        self.high_solubility_values = [0.00007, 0.00016, 0.1, 1000]
        self.medium_solubility_values = [0.00001, 0.00002, 0.00004, 0.00006]
        self.low_solubility_values = [0.000001, 0.000002, 0.0000025, 0.000005, 0.0000099]

    def test_all_solubility(self):
        self.assertTrue(all(SOLUBILITY['all'](v) for v in self.high_solubility_values))
        self.assertTrue(all(SOLUBILITY['all'](v) for v in self.medium_solubility_values))
        self.assertTrue(all(SOLUBILITY['all'](v) for v in self.low_solubility_values))

    def test_high_solubility(self):
        self.assertTrue(all(SOLUBILITY['high'](v) for v in self.high_solubility_values))
        self.assertFalse(any(SOLUBILITY['high'](v) for v in self.medium_solubility_values))
        self.assertFalse(any(SOLUBILITY['high'](v) for v in self.low_solubility_values))

    def test_medium_solubility(self):
        self.assertFalse(any(SOLUBILITY['medium'](v) for v in self.high_solubility_values))
        self.assertTrue(any(SOLUBILITY['medium'](v) for v in self.medium_solubility_values))
        self.assertFalse(any(SOLUBILITY['medium'](v) for v in self.low_solubility_values))

    def test_low_solubility(self):
        self.assertFalse(any(SOLUBILITY['low'](v) for v in self.high_solubility_values))
        self.assertFalse(any(SOLUBILITY['low'](v) for v in self.medium_solubility_values))
        self.assertTrue(all(SOLUBILITY['low'](v) for v in self.low_solubility_values))


class TanimotoUnitTest(unittest.TestCase):
    def setUp(self):
        self.fp1 = np.array([1, 2, 3], dtype=np.int32)
        self.fp2 = np.array([1, 2, 4], dtype=np.int32)
        self.fp3 = np.array([1, 2, 3, 4], dtype=np.int32)

    def match_test(self):
        self.assertEqual(tanimoto_coefficient(self.fp1, self.fp1), 1)
        self.assertEqual(tanimoto_distance(self.fp1, self.fp1), 0)

    def mismatch_test(self):
        self.assertGreater(tanimoto_coefficient(self.fp1, self.fp1), tanimoto_coefficient(self.fp1, self.fp2))
        self.assertGreater(tanimoto_coefficient(self.fp1, self.fp2), tanimoto_coefficient(self.fp1, self.fp3))
        self.assertGreater(tanimoto_coefficient(self.fp1, self.fp2), tanimoto_coefficient(self.fp2, self.fp3))

    def order_test(self):
        self.assertEqual(tanimoto_coefficient(self.fp1, self.fp2), tanimoto_coefficient(self.fp2, self.fp1))

    def different_size_test(self):
        self.assertEqual(tanimoto_coefficient(self.fp1, self.fp3), -1)
        self.assertEqual(tanimoto_distance(self.fp1, self.fp3), 2)

        self.assertEqual(tanimoto_coefficient(self.fp2, self.fp3), -1)
        self.assertEqual(tanimoto_distance(self.fp2, self.fp3), 2)

    def distance_is_one_minus_coefficient_test(self):
        self.assertAlmostEqual(tanimoto_distance(self.fp1, self.fp1), 1 - tanimoto_coefficient(self.fp1, self.fp1),
                               delta=0.00001)
        self.assertAlmostEqual(tanimoto_distance(self.fp1, self.fp2), 1 - tanimoto_coefficient(self.fp1, self.fp2),
                               delta=0.00001)
        self.assertAlmostEqual(tanimoto_distance(self.fp1, self.fp3), 1 - tanimoto_coefficient(self.fp1, self.fp3),
                               delta=0.00001)


class TestCaseWrapper:

    class ChemistryLibraryTestCaseWrapper(unittest.TestCase):

        def setUp(self):
            # Test molecule is L-Tryptophan (CHEBI:16828)
            self.inchi = "InChI=1S/C11H12N2O2/c12-9(11(14)15)5-7-6-13-10-4-2-1-3-8(7)10/h1-4,6,9,13H,5,12H2,(H,14,15)/t9-/m0/s1"
            self.inchi_key = "QIVBCDIJIAJPQS-VIFPVBQESA-N"
            self.smiles = "N[C@@H](Cc1c[nH]c2ccccc12)C(O)=O"
            self.formula = "C11H12N2O2"

        def _num_atoms(self, mol):
            raise NotImplementedError

        def _num_bonds(self, mol):
            raise NotImplementedError

        def _num_carbon(self, mol):
            raise NotImplementedError

        def _num_protons(self, mol):
            raise NotImplementedError

        def molecule_from_inchi_test(self):
            mol = self.library.inchi_to_molecule(self.inchi)
            self.assertEqual(self._num_atoms(mol), 27)
            self.assertEqual(self._num_bonds(mol), 28)
            self.assertEqual(self._num_carbon(mol), 11)
            self.assertEqual(self._num_protons(mol), 12)

        def inchi_to_inchi_key_test(self):
            inchi_key = self.library.inchi_to_inchi_key(self.inchi)
            self.assertEqual(self.inchi_key, inchi_key)


class OpenBabelTestCase(TestCaseWrapper.ChemistryLibraryTestCaseWrapper):
    def setUp(self):
        super(OpenBabelTestCase, self).setUp()
        self.library = openbabel

    def _num_atoms(self, mol):
        return len(mol.atoms)

    def _num_bonds(self, mol):
        return mol.OBMol.NumBonds()

    def _num_carbon(self, mol):
        return len([a for a in mol.atoms if a.atomicnum == CARBON_ATOMIC_NUMBER])

    def _num_protons(self, mol):
        return len([a for a in mol.atoms if a.atomicnum == HYDROGEN_ATOMIC_NUMBER])


class RdkitTestCase(TestCaseWrapper.ChemistryLibraryTestCaseWrapper):
    def setUp(self):
        super(RdkitTestCase, self).setUp()
        self.library = rdkit

    def _num_atoms(self, mol):
        return mol.GetNumAtoms()

    def _num_bonds(self, mol):
        return mol.GetNumBonds()

    def _num_carbon(self, mol):
        return len([a for a in mol.GetAtoms() if a.GetAtomicNum() == CARBON_ATOMIC_NUMBER])

    def _num_protons(self, mol):
        return len([a for a in mol.GetAtoms() if a.GetAtomicNum() == HYDROGEN_ATOMIC_NUMBER])
