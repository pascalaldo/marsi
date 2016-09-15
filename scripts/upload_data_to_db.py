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

from datetime import datetime


import os


from IProgress.progressbar import ProgressBar
from IProgress.widgets import Bar, ETA
from mongoengine import connect
from pandas import DataFrame

import openbabel as ob
from pybel import readfile

from marsi import processing
from marsi.io.data import kegg, drugbank, chebi, pubchem
from marsi.utils import data_dir, log_dir

ob.obErrorLog.SetOutputLevel(0)

summary = DataFrame(columns=["kegg", "chebi", "drugbank", "pubchem", "valid"])
progress = ProgressBar(maxval=len(kegg) + len(chebi) + len(drugbank) + len(pubchem),
                       widgets=["Adding metabolites: ", Bar(), ETA()])

connect(db="marsi-db")

progress.start()
i = 0


def fill_key(df, inchi_key, column, identifier):
    if inchi_key not in df.index.values:
        df.loc[inchi_key, column] = [identifier]
    elif isinstance(df.loc[inchi_key, column], float):
        df.loc[inchi_key, column] = [identifier]
    else:
        df.loc[inchi_key, column].append(identifier)

# import ChEBI data
for mol in readfile("sdf", os.path.join(data_dir, "chebi_lite_3star.sdf")):
    chebi_id = processing.mol_chebi_id(mol)
    if int(chebi_id.split(":")[1]) in chebi.compound_id.values:
        if not processing.has_radical(mol):
            inchi_key = processing.mol_to_inchi_key(mol)
            if len(inchi_key) > 0:
                processing.add_mol_to_db(mol, [chebi_id])
                fill_key(summary, inchi_key, 'chebi', chebi_id)
        i += 1
        progress.update(i)

# import DrugBank
for mol in readfile("sdf", os.path.join(data_dir, "drugbank_open_structures_5.0.1.sdf")):
    drugbank_id = processing.mol_drugbank_id(mol)
    if drugbank_id in drugbank.id.values:
        if not processing.has_radical(mol):
            inchi_key = processing.mol_to_inchi_key(mol)
            if len(inchi_key) > 0:
                processing.add_mol_to_db(mol, [drugbank_id])
                fill_key(summary, inchi_key, 'drugbank', drugbank_id)
        i += 1
        progress.update(i)


# import KEGG
for mol_file in os.listdir(os.path.join(data_dir, "kegg_mol_files")):
    if mol_file[-4:] == ".mol":
        kegg_id = mol_file[:-4]
        try:
            mol = next(readfile("mol", os.path.join(data_dir, "kegg_mol_files", mol_file)))
        except StopIteration:
            continue
        if not processing.has_radical(mol):
            inchi_key = processing.mol_to_inchi_key(mol)
            if len(inchi_key) > 0:
                processing.add_mol_to_db(mol, [kegg_id])
                fill_key(summary, inchi_key, 'kegg', kegg_id)
        i += 1
        progress.update(i)


# import PubChem
for mol_file in os.listdir(os.path.join(data_dir, "pubchem_sdf_files")):
    if mol_file[-4:] == ".sdf":
        pubchem_id = mol_file[:-4]
        try:
            mol = next(readfile("mol", os.path.join(data_dir, "pubchem_sdf_files", mol_file)))
        except StopIteration:
            continue
        if not processing.has_radical(mol):
            inchi_key = processing.mol_to_inchi_key(mol)
            if len(inchi_key) > 0:
                processing.add_mol_to_db(mol, [pubchem_id])
                fill_key(summary, inchi_key, 'pubchem', pubchem_id)
        i += 1
        progress.update(i)

progress.finish()

timestamp = datetime.now().strftime('%Y%m%d_%H_%M_%S')
summary.to_csv(os.path.join(log_dir, "%s_summary.csv" % timestamp))
