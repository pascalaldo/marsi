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

kegg = DataFrame(columns=['group', 'family', 'level', 'target', 'generic_name', 'name', 'drug_type', 'kegg_drug_id'])

with open(os.path.join(data_dir, "kegg_brite_08310.keg")) as kegg_data:
    group = None
    family = None
    generic_name = None
    level = None
    target = None
    i = 0
    for line in kegg_data:
        line = line.strip("\n")
        if line.startswith("A"):
            group = line[1:].strip("<b>").strip("<b/>")
        if group != "Enzymes":
            continue
            # if line.startswith("B"):
            #     family = line[1:].strip()
            # elif line.startswith("C"):
            #     level = line[1:].strip()
            # elif line.startswith("D"):
            #     target = line[1:].strip()
            # elif line.startswith("E"):
            #     generic_name = line[1:].strip()
            # elif line.startswith("F"):
            #     line = line[1:].strip()
            #     split = line.split()
            #     name = " ".join(split[1:-2])
            #     kegg.loc[i] = [group, family, level, target, generic_name, name, split[-1], split[0]]
            #     i += 1
        else:
            if line.startswith("B"):
                family = line[1:].strip()
                level = family
            elif line.startswith("C"):
                target = line[1:].strip()
            elif line.startswith("D"):
                generic_name = line[1:].strip()
            elif line.startswith("E"):
                line = line[1:].strip()
                split = line.split()
                name = " ".join(split[1:-2])
                kegg.loc[i] = [group, family, level, target, generic_name, name, split[-1], split[0]]
                i += 1
print("Found %i drugs acting on enzymes" % i)
kegg.to_csv(os.path.join(data_dir, "kegg_data.csv"))