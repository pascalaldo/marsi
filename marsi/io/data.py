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
from pandas import DataFrame

from marsi.utils import data_dir

__all__ = ["chebi", "drugbank"]

chebi = DataFrame.from_csv(os.path.join(data_dir, "chebi_analogues_filtered.csv"))

drugbank = DataFrame.from_csv(os.path.join(data_dir, "drugbank_open_vocabulary.csv"), sep=",", index_col=None)
drugbank.columns = ["id", "accessions", "common_names", "cas", "unii", "synonyms", "inchi_key"]

drugbank.fillna("", inplace=True)

drugbank['accessions'] = drugbank.accessions.apply(str.split, args=(" | ",))
drugbank['synonyms'] = drugbank.synonyms.apply(str.split, args=(" | ",))

pubchem = DataFrame.from_csv(os.path.join(data_dir, "pubchem_data.csv"))

kegg = DataFrame.from_csv(os.path.join(data_dir, "kegg_data.csv"))

# binding_db = read_csv(os.path.join(data_dir, "BindingDB_All.tsv"), sep="\t", error_bad_lines=False)
# with open(os.path.join(data_dir, 'binding_db_sane_columns.txt'), 'r') as columns_file:
#     line = next(columns_file)
#     binding_db.columns = line.split(", ")

# solubility = DataFrame.from_csv(os.path.join(data_dir, 'solubility.csv'))
# solubility.columns = ['log_measured', 'log_predicted', 'smiles']
