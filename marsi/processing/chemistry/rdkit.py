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
import math
from functools import lru_cache

try:
    import rdkit

    from rdkit import Chem
    from rdkit.Chem import MCS
except ImportError:
    class RdkitFail:
        def __dir__(self):
            return ["rdkit_not_available"]

        def __getattr__(self, item):
            return "RDKit is not installed"
    rdkit = RdkitFail()


@lru_cache(maxsize=256)
def inchi_to_molecule(inchi):
    """
    Returns a molecule from a InChI string.

    Arguments
    ---------
    inchi: str
        A valid string.

    Returns
    -------
    mol: rdkit.Chem.rdchem.Mol
        A molecule.
    """
    mol = Chem.MolFromInchi(inchi)
    mol = Chem.AddHs(mol)
    return mol


@lru_cache(maxsize=256)
def inchi_to_inchi_key(inchi):
    """
    Makes an InChI Key from a InChI string.

    Arguments
    ---------
    inchi: str
        A valid InChI string.

    Returns
    -------
    str
        A InChI key.
    """
    return Chem.InchiToInchiKey(inchi)


def mol_to_inchi_key(mol):
    """
    Makes an InChI Key from a pybel.Molecule.

    Arguments
    ---------
    mol: rdkit.Chem.rdchem.Mol
        A molecule.

    Returns
    -------
    str
        A InChI key.
    """
    return inchi_to_inchi_key(Chem.MolToInchi(mol))


def fingerprint(mol, fpformat='maccs'):
    """
    Returns the Fingerprint of the molecule.

    Arguments
    ---------
    mol: rdkit.Chem.rdchem.Mol
        A molecule.

    fpformat: str
        A valid fingerprint format.

    Returns
    -------
    Fingerprint
        A fingerprint
    """
    raise NotImplemented


def maximum_common_substructure(reference, molecule, match_rings=True, match_fraction=0.6, timeout=None):
    """
    Returns the Maximum Common Substructure (MCS) between two molecules.

    Arguments
    ---------
    reference: rdkit.Chem.Mol
        A molecule.
    molecule: rdkit.Chem.Mol
        Another molecule.
    match_rings: bool
        Force ring structure to match
    match_fraction: float
        Match is fraction of the reference atoms (default: 0.6)
    timeout: int

    Returns
    -------
    rdkit.Chem.MCS.MCSResult
    """

    assert isinstance(reference, rdkit.Chem.rdchem.Mol)

    min_num_atoms = math.ceil(reference.GetNumAtoms()) * match_fraction

    return MCS.FindMCS([reference, molecule], ringMatchesRingOnly=match_rings,
                       minNumAtoms=min_num_atoms, timeout=timeout, atomCompare='any')


def mcs_similarity(mcs_result, molecule, atoms_weight=0.5, bonds_weight=0.5):
    """
    Returns the Maximum Common Substructure (MCS) between two molecules.

    $$ atoms\_weight * (mcs_res.similar_atoms/mol.num_atoms) + bonds\_weight * (mcs_res.similar_bonds/mol.num_bonds) $$

    Arguments
    ---------
    mcs_result: rdkit.Chem.MCS.MCSResult
        The result of a Maximum Common Substructure run.
    molecule: rdkit.Chem.Mol
        A molecule.
    atoms_weight: float
        How much similar atoms matter.
    bonds_weight: float
        How much similar bonds matter.
    Returns
    -------
    float
    """

    assert isinstance(mcs_result, MCS.MCSResult)
    assert isinstance(molecule, Chem.rdchem.Mol)

    if mcs_result.completed:
        atoms_score = atoms_weight * (mcs_result.numAtoms / molecule.GetNumAtoms())
        bonds_score = bonds_weight * (mcs_result.numBonds / molecule.GetNumBonds())
        return atoms_score + bonds_score

    else:
        return .0


def structural_similarity(reference, molecule, atoms_weight=0.5, bonds_weight=0.5,
                          match_rings=True, match_fraction=0.6, timeout=None):
    """
    Returns a structural similarity based on the Maximum Common Substructure (MCS) between two molecules.

    $$ mcs_s(ref) * mcs_s(mol) $$

    Arguments
    ---------
    reference: rdkit.Chem.Mol
        The result of a Maximum Common Substructure run.
    molecule: rdkit.Chem.Mol
        A molecule.
    atoms_weight: float
        How much similar atoms matter.
    bonds_weight: float
        How much similar bonds matter.
    match_rings: bool
        Force ring structure to match
    match_fraction: float
        Match is fraction of the reference atoms (default: 0.6)
    timeout: int
    Returns
    -------
    float
    """

    mcs_res = maximum_common_substructure(reference, molecule, match_rings=match_rings,
                                          match_fraction=match_fraction, timeout=timeout)

    ref_similarity = mcs_similarity(mcs_res, reference, atoms_weight=atoms_weight, bonds_weight=bonds_weight)
    mol_similarity = mcs_similarity(mcs_res, molecule, atoms_weight=atoms_weight, bonds_weight=bonds_weight)
    return ref_similarity * mol_similarity
