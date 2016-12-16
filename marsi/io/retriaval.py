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
import bioservices
import os
from urllib.request import urlretrieve
import pubchempy as pbc

from IProgress import ProgressBar
from pandas import DataFrame

from marsi.utils import data_dir, gunzip
from ftplib import FTP
import requests

BIGG_BASE_URL = "http://bigg.ucsd.edu/static/namespace/"
DRUGBANK_BASE_URL = "https://www.drugbank.ca/"
CHEBI_FTP_URL = "ftp.ebi.ac.uk"
KEGG_BASE_URL = "http://www.genome.jp"
ZINC_BASE_URL = "http://zinc.docking.org/"


def retrieve_bigg_reactions(dest=os.path.join(data_dir, "bigg_models_reactions.txt")):
    """
    Retrieves bigg reactions file
    """
    bigg_reactions_file = "bigg_models_reactions.txt"
    urlretrieve(BIGG_BASE_URL+bigg_reactions_file, dest)

def retrieve_bigg_metabolites(dest=os.path.join(data_dir, "bigg_models_metabolites.txt")):
    """
    Retrieves bigg metabolites file
    """
    bigg_metabolites_file = "bigg_models_metabolites.txt"
    urlretrieve(BIGG_BASE_URL+bigg_metabolites_file, dest)


def retrieve_drugbank_open_structures(db_version="5.0.3", dest=os.path.join(data_dir, "drugbank_open_structures.sdf")):
    """
    Retrieves Drugbank Open Structures.

    Parameters
    ----------
    db_version: str
        The version of drugbank to retrieve
    """

    encoded_version = db_version.replace(".", "-")
    response = requests.get(DRUGBANK_BASE_URL + "releases/%s/downloads/all-open-structures" % encoded_version)
    with open(dest, 'w') as fhandle:
        fhandle.write(response.text)


def retrieve_drugbank_open_vocabulary(db_version="5.0.3", dest=os.path.join(data_dir,"drugbank_open_vocabulary.sdf")):
    """
    Retrieves Drugbank Open Vocabulary.

    Parameters
    ----------
    db_version: str
        The version of drugbank to retrieve
    """

    encoded_version = db_version.replace(".", "-")
    response = requests.get(DRUGBANK_BASE_URL + "releases/%s/downloads/all-drugbank-vocabulary" % encoded_version)
    with open(dest, 'w') as fhandle:
        fhandle.write(response.text)


def retrieve_chebi_structures(dest=os.path.join(data_dir, "chebi_lite_3star.sdf")):
    """
    Retrieves ChEBI sdf (lite version).
    """
    sdf_file = "ChEBI_lite_3star.sdf.gz"
    chebi_structures_file = dest + ".gz"
    ftp = FTP(CHEBI_FTP_URL)
    ftp.login()
    ftp.cwd('pub/databases/chebi/SDF')
    with open(chebi_structures_file, "wb") as structures_file:
        ftp.retrbinary("RETR %s" % sdf_file, structures_file.write)
    ftp.quit()
    gunzip(chebi_structures_file)


def retrieve_chebi_names(dest=os.path.join(data_dir, "chebi_names_3star.txt")):
    """
    Retrieves ChEBI names.
    """
    gz_file = "names_3star.tsv.gz"
    chebi_names_file = dest + ".gz"
    ftp = FTP(CHEBI_FTP_URL)
    ftp.login()
    ftp.cwd('pub/databases/chebi/Flat_file_tab_delimited')
    with open(chebi_names_file, "wb") as names_file:
        ftp.retrbinary("RETR %s" % gz_file, names_file.write)
    ftp.quit()
    gunzip(chebi_names_file)


def retrieve_chebi_relation(dest=os.path.join(data_dir, "chebi_relation_3star.tsv")):
    """
    Retrieves ChEBI relation data.
    """
    tsv_file = "relation_3star.tsv"
    ftp = FTP(CHEBI_FTP_URL)
    ftp.login()
    ftp.cwd('pub/databases/chebi/Flat_file_tab_delimited')
    with open(dest, "wb") as relation_file:
        ftp.retrbinary("RETR %s" % tsv_file, relation_file.write)
    ftp.quit()


def retrieve_chebi_vertice(dest=os.path.join(data_dir, "chebi_vertice_3star.tsv")):
    """
    Retrieves ChEBI vertice data.
    """
    tsv_file = "vertice_3star.tsv"
    ftp = FTP(CHEBI_FTP_URL)
    ftp.login()
    ftp.cwd('pub/databases/chebi/Flat_file_tab_delimited')
    with open(dest, "wb") as verice_file:
        ftp.retrbinary("RETR %s" % tsv_file, verice_file.write)
    ftp.quit()


def retrieve_kegg_brite(dest=os.path.join(data_dir, "kegg_brite_08310.keg")):
    """
    Retrieves KEGG Brite 08310 (Target-based Classification of Drugs)

    """
    urlretrieve(KEGG_BASE_URL + "/kegg-bin/download_htext?htext=br08310.keg&format=htext&filedir=", dest)


def retrieve_pubchem_mol_files(pubchem_ids, dest=data_dir):
    pubchem_files_path = os.path.join(dest, 'pubchem_sdf_files')

    if not os.path.isdir(pubchem_files_path):
        os.mkdir(pubchem_files_path)

    for pubchem_id in pubchem_ids:
        try:
            pbc.download('sdf', os.path.join(pubchem_files_path, '%i.sdf' % pubchem_id), int(pubchem_id))
        except IOError:
            # File already exists
            continue


def retrieve_kegg_mol_files(kegg, dest=data_dir):
    kegg_client = bioservices.kegg.KEGG()
    drug_ids = kegg.kegg_drug_id.unique()

    not_found = []

    kegg_mol_files_dir = os.path.join(dest, "kegg_mol_files")
    if not os.path.isdir(kegg_mol_files_dir):
        os.mkdir(kegg_mol_files_dir)

    for drug_id in drug_ids:
        file_name = os.path.join(kegg_mol_files_dir, "%s.mol" % drug_id)
        if not os.path.exists(file_name):
            with open(file_name, "w") as mol_file_handler:
                kegg_mol_data = kegg_client.get(drug_id, 'mol')
                if isinstance(kegg_mol_data, int) and kegg_mol_data == 404:
                    not_found.append(drug_id)
                elif len(kegg_mol_data.strip()) == 0:
                    not_found.append(drug_id)
                else:
                    mol_file_handler.write(kegg_mol_data)

    print("Not Found: %s" % (", ".join(not_found)))


def retrieve_zinc_properties(dest=os.path.join(data_dir, "zinc15_16_prop.tsv")):
    """
    Retrieves ZINC properties file:
    "All Clean
    As Subset #6, but without 'yuck' compounds"

    """
    urlretrieve(ZINC_BASE_URL + "db/bysubset/16/16_prop.xls", dest)
