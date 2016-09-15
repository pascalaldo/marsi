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
import re

from pandas import DataFrame
from marsi.utils import data_dir

pubchem = DataFrame(columns=["name", "molecular_weight", "formula", "uipac_name", "create_date", "compound_id"])

with open(os.path.join(data_dir, "pubchem_compound_analogs_antimetabolites.txt")) as pubchem_data:
    row = dict(name=None, molecular_weight=None, formula=None, uipac_name=None,
               create_date=None, compound_id=None)
    i = 0
    for line in pubchem_data:
        line = line.strip("\n")
        if len(line) == 0:
            if any(v for v in row.values()):
                pubchem.loc[i] = [row[k] for k in pubchem.columns]
                i += 1
            row = dict(name=None, molecular_weight=None, formula=None,
                       uipac_name=None, create_date=None, compound_id=None)
        elif re.match("^\d+\.*", line):
            row['name'] = line.split(". ", 1)[1].split("; ")[0]
        elif line.startswith("MW:"):
            match = re.match("MW:\s+(\d+\.\d+).+MF:\s(\w+)", line)
            row['molecular_weight'] = float(match.group(1))
            row['formula'] = match.group(2)
        elif line.startswith("IUPAC name:"):
            row['uipac_name'] = line[10:]
        elif line.startswith("Create Date:"):
            row['create_date'] = line[12:]
        elif line.startswith("CID:"):
            row['compound_id'] = int(line[5:])

pubchem['compound_id'] = pubchem.compound_id.apply(int)
pubchem.to_csv(os.path.join(data_dir, "pubchem_data.csv"))
