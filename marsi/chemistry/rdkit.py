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

import rdkit
from cachetools import cached, LRUCache

from rdkit import Chem
from rdkit.Chem import MCS, AllChem, MACCSkeys, EState

import time
import numpy as np
from bitarray import bitarray

from marsi.chemistry.common import monte_carlo_volume as mc_vol, inchi_key_lru_cache


lru_cache = LRUCache(maxsize=256)

periodic_table = Chem.GetPeriodicTable()


fps = ["maccs", "morgan2", "morgan3", "morgan4", "morgan5"]


@cached(lru_cache)
def inchi_to_molecule(inchi):
    """
    Returns a molecule from a InChI string.

    Parameters
    ----------
    inchi : str
        A valid string.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        A molecule.
    """
    mol = Chem.MolFromInchi(inchi)
    mol = Chem.AddHs(mol)
    return mol


def mol_to_molecule(file_or_molecule_desc, from_file=True):
    """
    Returns a molecule from a MOL file.

    Parameters
    ----------
    file_or_molecule_desc : str
        A valid MOL file path or a valid MOL string.
    from_file : bool
        If True tries to read the molecule from a file.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        A molecule.
    """
    if from_file:
        mol = Chem.MolFromMolFile(file_or_molecule_desc)
    else:
        mol = Chem.MolFromMolBlock(file_or_molecule_desc)
    mol = Chem.AddHs(mol)
    return mol


def sdf_to_molecule(file_or_molecule_desc, from_file=True):
    """
    Returns a molecule from a SDF file.

    Parameters
    ----------
    file_or_molecule_desc : str
        A valid sdf file path or a valid SDF string.
    from_file : bool
        If True tries to read the molecule from a file.

    Returns
    -------
    rdkit.Chem.rdchem.Mol
        A molecule.
    """

    if from_file:
        supplier = Chem.SDMolSupplier(file_or_molecule_desc)
    else:
        supplier = Chem.SDMolSupplier()
        supplier.SetData(file_or_molecule_desc, strictParsing=False)
    mol = next(supplier)
    mol = Chem.AddHs(mol)
    return mol


@cached(inchi_key_lru_cache)
def inchi_to_inchi_key(inchi):
    """
    Makes an InChI Key from a InChI string.

    Parameters
    ----------
    inchi : str
        A valid InChI string.

    Returns
    -------
    str
        A InChI key.
    """
    return Chem.InchiToInchiKey(inchi)


def mol_to_inchi_key(mol):
    """
    Makes an InChI Key from a Molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        A molecule.

    Returns
    -------
    str
        A InChI key.
    """
    return inchi_to_inchi_key(mol_to_inchi(mol))


def mol_to_inchi(mol):
    """
    Makes an InChI from a Molecule.

    Parameters
    ----------
    mol : rdkit.Chem.rdchem.Mol
        A molecule.

    Returns
    -------
    str
        A InChI.
    """
    return Chem.MolToInchi(mol)


def fingerprint(molecule, fpformat='maccs'):
    """
    Returns the Fingerprint of the molecule.

    Parameters
    ----------
    molecule : rdkit.Chem.rdchem.Mol
        A molecule.

    fpformat : str
        A valid fingerprint format.

    Returns
    -------
    Fingerprint
        rdkit.DataStructs.cDataStructs.ExplicitBitVect
    """
    assert isinstance(molecule, Chem.rdchem.Mol)
    assert isinstance(fpformat, str)
    fp = None

    if fpformat not in fps:
        raise AssertionError("'%s' is not a valid fingerprint format" % fpformat)
    if fpformat == "maccs":
        fp = MACCSkeys.FingerprintMol(molecule)
    elif fpformat == "morgan2":
        fp = AllChem.GetMorganFingerprintAsBitVect(molecule, 2)
    elif fpformat == "morgan3":
        fp = AllChem.GetMorganFingerprintAsBitVect(molecule, 3)
    elif fpformat == "morgan4":
        fp = AllChem.GetMorganFingerprintAsBitVect(molecule, 4)
    elif fpformat == "morgan5":
        fp = AllChem.GetMorganFingerprintAsBitVect(molecule, 5)

    return fp


def fingerprint_to_bits(fp, bits=1024):
    """
    Converts a pybel.Fingerprint into a binary array

    Parameters
    ----------
    fp : rdkit.DataStructs.cDataStructs.ExplicitBitVect
        A fingerprint molecule.
    bits : int
        Number of bits (default is 1024)
    Returns
    -------
    bitarray
    """
    bits_list = bitarray(bits)

    for i in range(fp.GenNumBits()):
        if fp.GetBit(i):
            bits_list[i-1] = 1

    return bits_list


def maximum_common_substructure(reference, molecule, match_rings=True, match_fraction=0.6, timeout=None):
    """
    Returns the Maximum Common Substructure (MCS) between two molecules.

    Parameters
    ----------
    reference : rdkit.Chem.Mol
        A molecule.
    molecule : rdkit.Chem.Mol
        Another molecule.
    match_rings : bool
        Force ring structure to match
    match_fraction : float
        Match is fraction of the reference atoms (default: 0.6)
    timeout: int

    Returns
    -------
    rdkit.Chem.MCS.MCSResult
        Maximum Common Substructure result.
    """

    assert isinstance(reference, rdkit.Chem.rdchem.Mol)

    min_num_atoms = math.ceil(reference.GetNumAtoms()) * match_fraction

    return MCS.FindMCS([reference, molecule], ringMatchesRingOnly=match_rings,
                       minNumAtoms=min_num_atoms, timeout=timeout, atomCompare='any', bondCompare="any")


def mcs_similarity(mcs_result, molecule, atoms_weight=0.5, bonds_weight=0.5):
    """
    Returns the Maximum Common Substructure (MCS) between two molecules.

    $$ atoms\_weight * (mcs_res.similar\_atoms/mol.num_atoms) + bonds\_weight *
    (mcs_res.similar\_bonds/mol.num\_bonds) $$

    Parameters
    ----------
    mcs_result : rdkit.Chem.MCS.MCSResult
        The result of a Maximum Common Substructure run.
    molecule : rdkit.Chem.Mol
        A molecule.
    atoms_weight : float
        How much similar atoms matter.
    bonds_weight : float
        How much similar bonds matter.
    Returns
    -------
    float
        Similarity value
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

    $$ mcs\_s(ref) * mcs\_s(mol) $$

    Parameters
    ----------
    reference : rdkit.Chem.Mol
        The result of a Maximum Common Substructure run.
    molecule : rdkit.Chem.Mol
        A molecule.
    atoms_weight : float
        How much similar atoms matter.
    bonds_weight : float
        How much similar bonds matter.
    match_rings : bool
        Force ring structure to match.
    match_fraction : float
        Match is fraction of the reference atoms (default: 0.6).
    timeout : int
    Returns
    -------
    float
        Similarity between reference and molecule.
    """

    mcs_res = maximum_common_substructure(reference, molecule, match_rings=match_rings,
                                          match_fraction=match_fraction, timeout=timeout)

    ref_similarity = mcs_similarity(mcs_res, reference, atoms_weight=atoms_weight, bonds_weight=bonds_weight)
    mol_similarity = mcs_similarity(mcs_res, molecule, atoms_weight=atoms_weight, bonds_weight=bonds_weight)
    return ref_similarity * mol_similarity


def monte_carlo_volume(molecule, coordinates=None, tolerance=1, max_iterations=10000, step_size=1000,
                       seed=time.time(), verbose=False, forcefield='mmff94', steps=100):
    """
    Adapted from:

    Simple Monte Carlo estimation of VdW molecular volume (in A^3)
    by Geoffrey Hutchison <geoffh@pitt.edu>

    https://github.com/ghutchis/hutchison-cluster

    Parameters
    ----------
    molecule : rdkit.Chem.rdchem.Mol
        A molecule from rdkit.
    coordinates : list
        A list of pre selected x,y,z coords. It must match the atoms order.
    tolerance : float
        The tolerance for convergence of the monte carlo abs(new_volume - volume) < tolerance
    max_iterations : int
        Number of iterations before the algorithm starts.
    step_size : int
        Number of points to add each step.
    seed : object
        A valid seed for random.
    verbose : bool
        Print debug information if True.
    forcefield : str
        The force field to get a 3D molecule. (only if it is not 3D already)
    steps : int
        The number of steps used for the force field to get a 3D molecule. (only if it is not 3D already)

    Returns
    -------
    float
        Molecule volume
    """

    assert isinstance(molecule, rdkit.Chem.rdchem.Mol)
    if coordinates is None:
        if len(molecule.GetConformers()) == 0:
            AllChem.EmbedMolecule(molecule)
            AllChem.UFFOptimizeMolecule(molecule)

        conformer = molecule.GetConformer(0)
        coordinates = []
        vdw_radii = []
        for index, atom in enumerate(molecule.GetAtoms()):
            pos = tuple(conformer.GetAtomPosition(index))
            radius = periodic_table.GetRvdw(atom.GetAtomicNum())
            coordinates.append(pos)
            vdw_radii.append(radius)
    coordinates = np.array(coordinates, dtype=np.float32)
    vdw_radii = np.array(vdw_radii, dtype=np.float32)
    return mc_vol(coordinates, vdw_radii, tolerance, max_iterations, step_size, seed, int(verbose))
