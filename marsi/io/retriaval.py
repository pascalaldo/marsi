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

from __future__ import absolute_import

import os
import zipfile
from ftplib import FTP
from io import BytesIO

import bioservices
import pubchempy as pbc
import requests
from IProgress import ProgressBar, Bar, ETA
from six.moves.urllib.request import urlretrieve

from marsi.utils import data_dir, gunzip

BIGG_BASE_URL = "http://bigg.ucsd.edu/static/namespace/"
DRUGBANK_BASE_URL = "https://www.drugbank.ca/"
CHEBI_FTP_URL = "ftp.ebi.ac.uk"
KEGG_BASE_URL = "http://www.genome.jp"
ZINC_BASE_URL = "http://zinc.docking.org/"

ZINC_STRUCTURES = ["16_p0.0.sdf.gz", "16_p0.1.sdf.gz", "16_p0.10.sdf.gz", "16_p0.100.sdf.gz", "16_p0.101.sdf.gz",
                   "16_p0.102.sdf.gz", "16_p0.103.sdf.gz", "16_p0.104.sdf.gz", "16_p0.105.sdf.gz", "16_p0.106.sdf.gz",
                   "16_p0.107.sdf.gz", "16_p0.108.sdf.gz", "16_p0.109.sdf.gz", "16_p0.11.sdf.gz", "16_p0.110.sdf.gz",
                   "16_p0.111.sdf.gz", "16_p0.112.sdf.gz", "16_p0.113.sdf.gz", "16_p0.114.sdf.gz", "16_p0.115.sdf.gz",
                   "16_p0.116.sdf.gz", "16_p0.117.sdf.gz", "16_p0.118.sdf.gz", "16_p0.119.sdf.gz", "16_p0.12.sdf.gz",
                   "16_p0.120.sdf.gz", "16_p0.121.sdf.gz", "16_p0.122.sdf.gz", "16_p0.123.sdf.gz", "16_p0.124.sdf.gz",
                   "16_p0.125.sdf.gz", "16_p0.126.sdf.gz", "16_p0.13.sdf.gz", "16_p0.14.sdf.gz", "16_p0.15.sdf.gz",
                   "16_p0.16.sdf.gz", "16_p0.17.sdf.gz", "16_p0.18.sdf.gz", "16_p0.19.sdf.gz", "16_p0.2.sdf.gz",
                   "16_p0.20.sdf.gz", "16_p0.21.sdf.gz", "16_p0.22.sdf.gz", "16_p0.23.sdf.gz", "16_p0.24.sdf.gz",
                   "16_p0.25.sdf.gz", "16_p0.26.sdf.gz", "16_p0.27.sdf.gz", "16_p0.28.sdf.gz", "16_p0.29.sdf.gz",
                   "16_p0.3.sdf.gz", "16_p0.30.sdf.gz", "16_p0.31.sdf.gz", "16_p0.32.sdf.gz", "16_p0.33.sdf.gz",
                   "16_p0.34.sdf.gz", "16_p0.35.sdf.gz", "16_p0.36.sdf.gz", "16_p0.37.sdf.gz", "16_p0.38.sdf.gz",
                   "16_p0.39.sdf.gz", "16_p0.4.sdf.gz", "16_p0.40.sdf.gz", "16_p0.41.sdf.gz", "16_p0.42.sdf.gz",
                   "16_p0.43.sdf.gz", "16_p0.44.sdf.gz", "16_p0.45.sdf.gz", "16_p0.46.sdf.gz", "16_p0.47.sdf.gz",
                   "16_p0.48.sdf.gz", "16_p0.49.sdf.gz", "16_p0.5.sdf.gz", "16_p0.50.sdf.gz", "16_p0.51.sdf.gz",
                   "16_p0.52.sdf.gz", "16_p0.53.sdf.gz", "16_p0.54.sdf.gz", "16_p0.55.sdf.gz", "16_p0.56.sdf.gz",
                   "16_p0.57.sdf.gz", "16_p0.58.sdf.gz", "16_p0.59.sdf.gz", "16_p0.6.sdf.gz", "16_p0.60.sdf.gz",
                   "16_p0.61.sdf.gz", "16_p0.62.sdf.gz", "16_p0.63.sdf.gz", "16_p0.64.sdf.gz", "16_p0.65.sdf.gz",
                   "16_p0.66.sdf.gz", "16_p0.67.sdf.gz", "16_p0.68.sdf.gz", "16_p0.69.sdf.gz", "16_p0.7.sdf.gz",
                   "16_p0.70.sdf.gz", "16_p0.71.sdf.gz", "16_p0.72.sdf.gz", "16_p0.73.sdf.gz", "16_p0.74.sdf.gz",
                   "16_p0.75.sdf.gz", "16_p0.76.sdf.gz", "16_p0.77.sdf.gz", "16_p0.78.sdf.gz", "16_p0.79.sdf.gz",
                   "16_p0.8.sdf.gz", "16_p0.80.sdf.gz", "16_p0.81.sdf.gz", "16_p0.82.sdf.gz", "16_p0.83.sdf.gz",
                   "16_p0.84.sdf.gz", "16_p0.85.sdf.gz", "16_p0.86.sdf.gz", "16_p0.87.sdf.gz", "16_p0.88.sdf.gz",
                   "16_p0.89.sdf.gz", "16_p0.9.sdf.gz", "16_p0.90.sdf.gz", "16_p0.91.sdf.gz", "16_p0.92.sdf.gz",
                   "16_p0.93.sdf.gz", "16_p0.94.sdf.gz", "16_p0.95.sdf.gz", "16_p0.96.sdf.gz", "16_p0.97.sdf.gz",
                   "16_p0.98.sdf.gz", "16_p0.99.sdf.gz", "16_p1.0.sdf.gz", "16_p1.1.sdf.gz", "16_p1.10.sdf.gz",
                   "16_p1.11.sdf.gz", "16_p1.12.sdf.gz", "16_p1.13.sdf.gz", "16_p1.14.sdf.gz", "16_p1.15.sdf.gz",
                   "16_p1.16.sdf.gz", "16_p1.17.sdf.gz", "16_p1.18.sdf.gz", "16_p1.19.sdf.gz", "16_p1.2.sdf.gz",
                   "16_p1.20.sdf.gz", "16_p1.21.sdf.gz", "16_p1.22.sdf.gz", "16_p1.3.sdf.gz", "16_p1.4.sdf.gz",
                   "16_p1.5.sdf.gz", "16_p1.6.sdf.gz", "16_p1.7.sdf.gz", "16_p1.8.sdf.gz", "16_p1.9.sdf.gz"]

ZINC_SUBSET_16_BASE = "http://zinc.docking.org/db/bysubset/16"


def retrieve_bigg_reactions(dest=os.path.join(data_dir, "bigg_models_reactions.txt")):
    """
    Retrieves bigg reactions file
    """
    bigg_reactions_file = "bigg_models_reactions.txt"
    urlretrieve(BIGG_BASE_URL + bigg_reactions_file, dest)


def retrieve_bigg_metabolites(dest=os.path.join(data_dir, "bigg_models_metabolites.txt")):
    """
    Retrieves bigg metabolites file
    """
    bigg_metabolites_file = "bigg_models_metabolites.txt"
    urlretrieve(BIGG_BASE_URL + bigg_metabolites_file, dest)


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
    with zipfile.ZipFile(BytesIO(response.content)) as zip_ref:
        zip_ref.extractall(data_dir)

    os.rename(os.path.join(data_dir, "open structures.sdf"), dest)


def retrieve_drugbank_open_vocabulary(db_version="5.0.3", dest=os.path.join(data_dir, "drugbank_open_vocabulary.csv")):
    """
    Retrieves Drugbank Open Vocabulary.

    Parameters
    ----------
    db_version: str
        The version of drugbank to retrieve
    """

    encoded_version = db_version.replace(".", "-")
    response = requests.get(DRUGBANK_BASE_URL + "releases/%s/downloads/all-drugbank-vocabulary" % encoded_version)
    with zipfile.ZipFile(BytesIO(response.content)) as zip_ref:
        zip_ref.extractall(data_dir)

    os.rename(os.path.join(data_dir, "drugbank vocabulary.csv"), dest)


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
    """
    Retrieves SDF Files from PubChem.
    """
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
    """
    Retrieves KEGG MOL Files using KEGG REST API.
    """
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


def retrieve_zinc_properties(dest=os.path.join(data_dir, "zinc_16_prop.tsv")):
    """
    Retrieves ZINC properties file:
    "All Clean
    As Subset #6, but without 'yuck' compounds"

    """
    urlretrieve(ZINC_BASE_URL + "db/bysubset/16/16_prop.xls", dest)


def retrieve_zinc_structures(dest=os.path.join(data_dir, "zinc_16.sdf.gz")):
    """
    Retrieves ZINC structures file:
    "All Clean
    As Subset #6, but without 'yuck' compounds"

    """
    pbar = ProgressBar(maxval=len(ZINC_STRUCTURES), widgets=["Downloading Zinc Structures 16: ", Bar(), ETA()])
    with open(dest, 'wb') as output_file:
        for sdf_file in pbar(ZINC_STRUCTURES):
            response = requests.get(ZINC_SUBSET_16_BASE + "/" + sdf_file)
            output_file.write(response.content)
