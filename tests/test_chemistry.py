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

import numpy as np
import pytest

from marsi.chemistry import openbabel, rdkit
from marsi.chemistry.common import SOLUBILITY, tanimoto_coefficient, tanimoto_distance
from marsi.chemistry.molecule import Molecule

TEST_DIR = os.path.dirname(__file__)

MOL_VOLUMES = {
    "Diphenylketene": 173.769,
    "Acetate": None,
    "Cobamamide": None
}

MOL_RINGS = {
    "Diphenylketene": 2,
    "Acetate": 0,
    "Cobamamide": 15
}

MOL_BONDS = {
    "Diphenylketene": 26,
    "Acetate": 6,
    "Cobamamide": 223
}

MOL_ATOMS = {
    "Diphenylketene": 25,
    "Acetate": 7,
    "Cobamamide": 209
}

MOL_CARBONS = {
    "Diphenylketene": 14,
    "Acetate": 2,
    "Cobamamide": 72
}

molecules = list(MOL_CARBONS.keys())


CARBON_ATOMIC_NUMBER = 6
HYDROGEN_ATOMIC_NUMBER = 1

INCHI = "InChI=1S/C11H12N2O2/c12-9(11(14)15)5-7-6-13-10-4-2-1-3-8(7)10/h1-4,6,9,13H,5,12H2,(H,14,15)/t9-/m0/s1"
INCHI_KEY = "QIVBCDIJIAJPQS-VIFPVBQESA-N"


class openbabel_handler(object):
    @staticmethod
    def num_atoms(mol):
        return len(mol.atoms)

    @staticmethod
    def num_bonds(mol):
        return mol.OBMol.NumBonds()

    @staticmethod
    def num_carbon(mol):
        return len([a for a in mol.atoms if a.atomicnum == CARBON_ATOMIC_NUMBER])

    @staticmethod
    def num_protons(mol):
        return len([a for a in mol.atoms if a.atomicnum == HYDROGEN_ATOMIC_NUMBER])


class rdkit_handler(object):
    @staticmethod
    def num_atoms(mol):
        return mol.GetNumAtoms()

    @staticmethod
    def num_bonds(mol):
        return mol.GetNumBonds()

    @staticmethod
    def num_carbon(mol):
        return len([a for a in mol.GetAtoms() if a.GetAtomicNum() == CARBON_ATOMIC_NUMBER])

    @staticmethod
    def num_protons(mol):
        return len([a for a in mol.GetAtoms() if a.GetAtomicNum() == HYDROGEN_ATOMIC_NUMBER])


@pytest.fixture(params=['rdkit', 'openbabel'])
def chemlib(request):
    if request.param == 'rdkit':
        return rdkit, rdkit_handler
    elif request.param == 'openbabel':
        return openbabel, openbabel_handler
    else:
        raise ValueError("Invalid param %s" % request.param)


def test_solubility_thresholds():
    high_solubility_values = [0.00007, 0.00016, 0.1, 1000]
    medium_solubility_values = [0.00001, 0.00002, 0.00004, 0.00006]
    low_solubility_values = [0.000001, 0.000002, 0.0000025, 0.000005, 0.0000099]

    assert all(SOLUBILITY['all'](v) for v in high_solubility_values)
    assert all(SOLUBILITY['all'](v) for v in medium_solubility_values)
    assert all(SOLUBILITY['all'](v) for v in low_solubility_values)

    assert all(SOLUBILITY['high'](v) for v in high_solubility_values)
    assert not any(SOLUBILITY['high'](v) for v in medium_solubility_values)
    assert not any(SOLUBILITY['high'](v) for v in low_solubility_values)

    assert not any(SOLUBILITY['medium'](v) for v in high_solubility_values)
    assert all(SOLUBILITY['medium'](v) for v in medium_solubility_values)
    assert not any(SOLUBILITY['medium'](v) for v in low_solubility_values)

    assert not any(SOLUBILITY['low'](v) for v in high_solubility_values)
    assert not any(SOLUBILITY['low'](v) for v in medium_solubility_values)
    assert all(SOLUBILITY['low'](v) for v in low_solubility_values)


def test_tanimoto_coefficient(benchmark):
    fp1 = np.array([1, 2, 3], dtype=np.int32)
    fp2 = np.array([1, 2, 4], dtype=np.int32)
    fp3 = np.array([1, 2, 3, 4], dtype=np.int32)
    coefficient = benchmark(tanimoto_coefficient, fp1, fp1)
    assert coefficient == 1

    assert tanimoto_coefficient(fp1, fp3) == -1
    assert tanimoto_coefficient(fp2, fp3) == -1

    assert tanimoto_coefficient(fp1, fp1) > tanimoto_coefficient(fp1, fp2)
    assert tanimoto_coefficient(fp1, fp2) > tanimoto_coefficient(fp1, fp3)
    assert tanimoto_coefficient(fp1, fp1) > tanimoto_coefficient(fp2, fp3)

    assert tanimoto_coefficient(fp1, fp2) == tanimoto_coefficient(fp2, fp1)


def test_tanimoto_distance(benchmark):
    fp1 = np.array([1, 2, 3], dtype=np.int32)
    fp2 = np.array([1, 2, 4], dtype=np.int32)
    fp3 = np.array([1, 2, 3, 4], dtype=np.int32)
    coefficient = benchmark(tanimoto_distance, fp1, fp1)
    assert coefficient == 0

    assert tanimoto_distance(fp1, fp3) == 2
    assert tanimoto_distance(fp2, fp3) == 2


def test_tanimoto_distance_vs_coefficient():
    fp1 = np.array([1, 2, 3], dtype=np.int32)
    fp2 = np.array([1, 2, 4], dtype=np.int32)
    fp3 = np.array([1, 2, 3, 4], dtype=np.int32)
    assert tanimoto_distance(fp1, fp1) == pytest.approx(1 - tanimoto_coefficient(fp1, fp1), 1e-6)
    assert tanimoto_distance(fp1, fp2) == pytest.approx(1 - tanimoto_coefficient(fp1, fp2), 1e-6)
    assert tanimoto_distance(fp1, fp3) == pytest.approx(1 - tanimoto_coefficient(fp1, fp3), 1e-6)


def test_molecule_from_inchi_test(chemlib, benchmark):
    mol = benchmark(chemlib[0].inchi_to_molecule, INCHI)
    assert chemlib[1].num_atoms(mol) == 27
    assert chemlib[1].num_bonds(mol) == 28
    assert chemlib[1].num_carbon(mol) == 11
    assert chemlib[1].num_protons(mol) == 12


def test_inchi_to_inchi_key(chemlib, benchmark):
    inchi_key = benchmark(chemlib[0].inchi_to_inchi_key, INCHI)
    assert inchi_key == INCHI_KEY


@pytest.fixture(params=['large', 'medium', 'small'])
def inchi(request):
    if request.param == "large":
        # Raspberry ellagitannin
        return "InChI=1S/C116H76O74/c117-30-1-18(2-31(118)61(30)130)100(158)189-115-98-95(183-107(165)" \
               "25-10-38(125)66(135)76(145)50(25)53-28(110(168)186-98)13-41(128)69(138)79(53)148)90-45" \
               "(177-115)16-173-103(161)21-6-34(121)71(140)81(150)55(21)57-59(112(170)180-90)92(87(156)" \
               "85(154)83(57)152)175-43-4-19(3-32(119)62(43)131)101(159)190-116-99-96(184-108(166)26-11" \
               "-39(126)67(136)77(146)51(26)54-29(111(169)187-99)14-42(129)70(139)80(54)149)91-46(178-116)" \
               "17-174-104(162)22-7-35(122)72(141)82(151)56(22)58-60(113(171)181-91)93(88(157)86(155)84(58)" \
               "153)188-114-97-94(182-106(164)24-9-37(124)65(134)75(144)49(24)52-27(109(167)185-97)12-40(127)" \
               "68(137)78(52)147)89-44(176-114)15-172-102(160)20-5-33(120)63(132)73(142)47(20)48-23(105(163)" \
               "179-89)8-36(123)64(133)74(48)143/h1-14,44-46,89-91,94-99,114-157H,15-17H2/t44-,45-,46-,89-,90-" \
               ",91-,94+,95+,96+,97-,98-,99-,114+,115+,116+/m1/s1"
    elif request.param == "medium":
        # NADP
        return "InChI=1S/C21H28N7O17P3/c22-17-12-19(25-7-24-17)28(8-26-12)21-16(44-46(33,34)35)14(30)11(43-21)" \
               "6-41-48(38,39)45-47(36,37)40-5-10-13(29)15(31)20(42-10)27-3-1-2-9(4-27)18(23)32/h1-4,7-8,10-11" \
               ",13-16,20-21,29-31H,5-6H2,(H7-,22,23,24,25,32,33,34,35,36,37,38,39)/p+1/t10-,11-,13-,14-,15-,16" \
               "-,20-,21-/m1/s1"

    elif request.param == "small":
        # Pyruvate
        return "InChI=1S/C3H4O3/c1-2(4)3(5)6/h1H3,(H,5,6)/p-1"


def test_structural_similarity_is_1(inchi, benchmark):
    mol = rdkit.inchi_to_molecule(inchi)
    similarity = benchmark(rdkit.structural_similarity, mol, mol)
    assert similarity == 1


def test_structural_similarity_to_glucose(inchi, benchmark):
    mol = rdkit.inchi_to_molecule(inchi)
    glucose = rdkit.inchi_to_molecule("InChI=1S/C6H12O6/c7-1-2-3(8)4(9)5(10)6(11)12-2/h2-11H,1H2/t2-,3-,4+,5-,6-/m1/s1")
    similarity = benchmark(rdkit.structural_similarity, glucose, mol)
    assert similarity < 1


def test_mol_to_inchi(chemlib, inchi, benchmark):
    mol = benchmark(chemlib[0].inchi_to_molecule, inchi)
    inchi_ = chemlib[0].mol_to_inchi(mol)
    assert inchi_ == inchi
    assert chemlib[0].inchi_to_inchi_key(inchi_) == chemlib[0].mol_to_inchi_key(mol)


class Mol3D(object):
    def __init__(self, molecule, volume):
        self.molecule = molecule
        self.volume = volume


@pytest.fixture(params=molecules)
def mol3d(request):
    molecule = os.path.join(TEST_DIR, "fixtures", "%s.sdf" % request.param)
    return Mol3D(molecule, MOL_VOLUMES[request.param])


@pytest.fixture(params=molecules)
def molecule(request):
    mol_path = os.path.join(TEST_DIR, "fixtures", "%s.sdf" % request.param)
    ob_mol = openbabel.sdf_to_molecule(mol_path)
    rd_mol = rdkit.sdf_to_molecule(mol_path)
    mol = Molecule(ob_mol, rd_mol)
    setattr(mol, 'id', request.param)
    return mol


@pytest.mark.skip("Not working as expected - also not part of the main workflow")
def test_volume(chemlib, mol3d, benchmark):
    molecule = chemlib[0].sdf_to_molecule(mol3d.molecule)
    volume = benchmark(chemlib[0].monte_carlo_volume, molecule, max_iterations=1000)
    assert mol3d.volume * .9 <= volume <= mol3d.volume * 1.1


def test_atom_count(molecule):
    assert molecule.num_atoms == MOL_ATOMS[molecule.id]


def test_bond_count(molecule):
    assert molecule.num_bonds == MOL_BONDS[molecule.id]


def test_ring_count(molecule):
    assert molecule.num_rings == MOL_RINGS[molecule.id]
