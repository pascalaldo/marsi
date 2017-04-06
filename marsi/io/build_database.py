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
import pybel
from pybel import readfile

from numpy import nan

from marsi.io.mongodb import Metabolite
from marsi.io.mongodb import Reference
from marsi.chemistry import openbabel


def build_database(data, data_dir):
    """
    Builds then Molecules database.
    It requires that the input files have been downloaded.

    See Also
    --------
    Initialization documentation.

    Parameters
    ----------
    data : module
        The marsi.io.data module
    data_dir : str
        The path to where data is stored.

    """
    i = 0
    chebi_structures_file = os.path.join(data_dir, "chebi_lite_3star.sdf")
    i = upload_chebi_entries(chebi_structures_file, data.chebi, i=i)
    print("Added %i" % i)
    drugbank_structures_file = os.path.join(data_dir, "drugbank_open_structures.sdf")
    i = upload_drugbank_entries(drugbank_structures_file, data.drugbank, i=i)
    print("Added %i" % i)
    kegg_mol_files_dir = os.path.join(data_dir, "kegg_mol_files")
    i = upload_kegg_entries(kegg_mol_files_dir, data.kegg, i=i)
    print("Added %i" % i)
    pubchem_sdf_files_dir = os.path.join(data_dir, "pubchem_sdf_files")
    i = upload_pubchem_entries(pubchem_sdf_files_dir, data.pubchem, i=i)
    print("Added %i" % i)
    zinc_data_file = os.path.join(data_dir, "zinc_16.sdf.gz")
    i = upload_zin_entries(zinc_data_file, i=i)
    print("Added %i" % i)
    return i


def _add_molecule(mol, synonyms, database, identifier, is_analog):
    """
    Add a molecule to the database. It checks for radicals, and it only adds complete molecules.

    Parameters
    ----------
    mol : pybel.Molecule
        A molecule parsed by openbabel.
    synonyms : list
        A list of strings with common names for the molecule.
    database : str
        The database from where this molecule was retrieved.
    identifier : str
        The molecule identifier at database.
    is_analog : bool
        If the metabolite was labled as an analog.

    """
    if not openbabel.has_radical(mol):
        inchi_key = openbabel.mol_to_inchi_key(mol)
        if len(inchi_key) > 0:
            reference = Reference.add_reference(database, identifier)
            Metabolite.from_molecule(mol, [reference], synonyms, is_analog)


def upload_chebi_entries(chebi_structures_file, chebi_data, i=0):
    """
    Import ChEBI data
    """
    for mol in readfile("sdf", chebi_structures_file):
        chebi_id = openbabel.mol_chebi_id(mol)
        chebi_id_int = int(chebi_id.split(":")[1])
        chebi_rows = chebi_data.query('compound_id == @chebi_id_int')
        if len(chebi_rows) > 0:
            synonyms = list(chebi_rows.name)
            _add_molecule(mol, synonyms, 'chebi', chebi_id, True)

        i += 1
    return i


def upload_drugbank_entries(drugbank_structures_file, drugbank_data, i=0):
    """
    Import DrugBank
    """
    for mol in readfile("sdf", drugbank_structures_file):
        drugbank_id = openbabel.mol_drugbank_id(mol)
        drugbank_rows = drugbank_data.query("id == @drugbank_id")
        if len(drugbank_rows) > 0:
            _add_molecule(mol, drugbank_rows.iloc[0].synonyms[0], 'drugbank', drugbank_id, False)

        i += 1
    return i


def upload_kegg_entries(kegg_mol_files_dir, kegg_data, i=0):
    """
    Import KEGG
    """
    for mol_file in os.listdir(kegg_mol_files_dir):
        if mol_file[-4:] == ".mol":
            kegg_id = mol_file[:-4]
            try:
                mol = next(readfile("mol", os.path.join(kegg_mol_files_dir, mol_file)))
            except StopIteration:
                continue
            rows = kegg_data.query("kegg_drug_id == @kegg_id")
            synonyms = set(rows.generic_name.values.tolist() + rows.name.values.tolist())
            if None in synonyms:
                synonyms.remove(None)
            if nan in synonyms:
                synonyms.remove(nan)
            try:
                _add_molecule(mol, synonyms, 'kegg', kegg_id, False)
            except Exception as e:
                print(synonyms)
                raise e

        i += 1
    return i


def upload_pubchem_entries(pubchem_sdf_files_dir, pubchem_data, i=0):
    """
    Import PubChem
    """

    for sdf_file in os.listdir(pubchem_sdf_files_dir):
        if sdf_file[-4:] == ".sdf":
            pubchem_id = sdf_file[:-4]
            try:
                mol = next(readfile("mol", os.path.join(pubchem_sdf_files_dir, sdf_file)))
            except StopIteration:
                continue
            rows = pubchem_data.query("compound_id == @pubchem_id")
            synonyms = set(rows.name.values.tolist() + rows.uipac_name.values.tolist())
            if None in synonyms:
                synonyms.remove(None)

            _add_molecule(mol, synonyms, 'pubchem', pubchem_id, True)

        i += 1

    return i


def upload_zin_entries(zinc_data_file, i=0):
    """
    Add ZINC
    """
    if os.path.isfile(zinc_data_file):
        zinc = pybel.readfile('sdf', zinc_data_file)
        for molecule in zinc:
            _add_molecule(molecule, [], 'zinc', molecule.title, False)
            i += 1

        return i
