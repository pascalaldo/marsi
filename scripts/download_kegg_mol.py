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
import bioservices
from IProgress.progressbar import ProgressBar
from IProgress.widgets import Bar, Percentage, ETA

from marsi.io.data import kegg
from marsi.utils import data_dir

kegg_client = bioservices.kegg.KEGG()

drug_ids = kegg.kegg_drug_id.unique()
progress = ProgressBar(widgets=["Downloading: ", Bar(), Percentage(), " ", ETA()], maxval=len(drug_ids))

not_found = []

for drug_id in progress(drug_ids):
    file_name = os.path.join(data_dir, "kegg_mol_files", "%s.mol" % drug_id)
    if not os.path.exists(file_name):
        with open(file_name, "w") as mol_file_handler:
            kegg_mol_file = kegg_client.get(drug_id, 'mol')
            if isinstance(kegg_mol_file, int) and kegg_mol_file == 404:
                not_found.append(drug_id)
            elif len(kegg_mol_file.strip()) == 0:
                not_found.append(drug_id)
            else:
                mol_file_handler.write(kegg_mol_file)

print("Not Found: %s" % (", ".join(not_found)))
