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

import pybel

from marsi.io.db import Metabolite


def add_mol_to_db(mol, references):
    try:
        metabolite = Metabolite.get(inchi_key=mol_to_inchi_key(mol))
        for reference in references:
            if reference not in metabolite.references:
                metabolite.references.append(reference)
    except KeyError:
        metabolite = Metabolite(inchi_key=mol_to_inchi_key(mol),
                                inchi=mol_to_inchi(mol),
                                references=references)
    metabolite.save()

    return metabolite


def has_radical(mol):
    """
    Finds if a pybel.Molecule has Radicals.
    Radicals have an atomic number of 0.

    Arguments
    ---------
    mol: pybel.Molecule
        A molecule.

    Returns
    -------
    bool
        True if there are any radicals.
    """
    return any(a.atomicnum == 0 for a in mol.atoms)


def mol_to_inchi(mol):
    """
    Makes an InChI from a pybel.Molecule.

    Arguments
    ---------
    mol: pybel.Molecule
        A molecule.

    Returns
    -------
    str
        A InChI string.
    """
    return mol.write(format="inchi", opt=dict(errorlevel=0)).strip()


def mol_to_inchi_key(mol):
    """
    Makes an InChI Key from a pybel.Molecule.

    Arguments
    ---------
    mol: pybel.Molecule
        A molecule.

    Returns
    -------
    str
        A InChI key.
    """
    return mol.write(format="inchikey", opt=dict(errorlevel=0)).strip()


def mol_drugbank_id(mol):
    """
    Returns the DrugBank ID from the molecule data.

    Arguments
    ---------
    mol: pybel.Molecule
        A molecule.

    Returns
    -------
    str
       DrugBank ID
    """
    return mol.data['DRUGBANK_ID'].strip()


def mol_pubchem_id(mol):
    """
    Returns the PubChem Compound ID from the molecule data.

    Arguments
    ---------
    mol: pybel.Molecule
        A molecule.

    Returns
    -------
    str
       PubChem Compound ID
    """
    return mol.data['PUBCHEM_COMPOUND_CID'].strip()


def mol_chebi_id(mol):
    """
    Returns the ChEBI ID from the molecule data.

    Arguments
    ---------
    mol: pybel.Molecule
        A molecule.

    Returns
    -------
    str
       ChEBI ID
    """
    return mol.data['ChEBI ID'].strip()


def fingerprint(mol, fpformat='maccs'):
    """
    Returns the ChEBI ID from the molecule data.

    Arguments
    ---------
    mol: pybel.Molecule
        A molecule.

    fpformat: str
        A valid fingerprint format (see pybel.fps)

    Returns
    -------
    pybel.Fingerprint
        A fingerprint
    """
    if fpformat not in pybel.fps:
        raise AssertionError("'%s' is not a valid fingerprint format" % fpformat)
    mol.calcfp(fptype=fpformat)
