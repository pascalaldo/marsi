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
import re
from pandas import DataFrame


def parse_kegg_brite(brite_file):
    kegg = DataFrame(columns=['group', 'family', 'level', 'target', 'generic_name',
                              'name', 'drug_type', 'kegg_drug_id'])

    with open(brite_file) as kegg_data:
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
    return kegg


def parse_chebi_data(chebi_names_file, chebi_vertice_file, chebi_relation_file):
    chebi_names = DataFrame.from_csv(chebi_names_file, sep="\t")
    chebi_names.fillna("", inplace=True)
    chebi_names.index.name = "id"

    chebi_names.columns = map(str.lower, chebi_names.columns)
    chebi_names.drop_duplicates('compound_id', keep='last', inplace=True)
    chebi_names['adapted'] = chebi_names.adapted.apply(lambda v: v == "T")

    chebi_analogues = chebi_names[chebi_names.name.str.contains('analog')]
    chebi_antimetabolite = chebi_names[chebi_names.compound_id == 35221]

    chebi_relations = DataFrame.from_csv(chebi_relation_file, sep="\t")
    chebi_relations.columns = map(str.lower, chebi_relations.columns)
    chebi_relations.index.name = "id"

    chebi_vertices = DataFrame.from_csv(chebi_vertice_file, sep="\t")
    chebi_vertices.columns = map(str.lower, chebi_vertices.columns)
    chebi_vertices.index.name = "id"

    def retrieve_child_id(compound_id):
        return chebi_vertices.loc[compound_id, 'compound_child_id']

    chebi_relations['init_compound_id'] = chebi_relations.init_id.apply(retrieve_child_id)
    chebi_relations['final_compound_id'] = chebi_relations.final_id.apply(retrieve_child_id)

    chebi_is_a = chebi_relations[chebi_relations['type'] == 'is_a']
    chebi_has_role = chebi_relations[chebi_relations['type'] == 'has_role']

    def recursive_search(roots, relations, universe, aggregated, forward=True):
        aggregated = aggregated.append(roots, ignore_index=True)
        if forward:
            filtered = relations[relations.init_compound_id.isin(roots.compound_id)]
            roots = universe[universe.compound_id.isin(filtered.final_compound_id)]
        else:
            filtered = relations[relations.final_compound_id.isin(roots.compound_id)]
            roots = universe[universe.compound_id.isin(filtered.init_compound_id)]

        if len(roots) > 0:
            aggregated, roots = recursive_search(roots, relations, universe, aggregated, forward)

        return aggregated, roots

    data = DataFrame(columns=chebi_names.columns)
    anti = DataFrame(columns=chebi_names.columns)

    data, _ = recursive_search(chebi_analogues, chebi_is_a, chebi_names, data, True)
    data, _ = recursive_search(chebi_antimetabolite, chebi_is_a, chebi_names, data, True)

    anti, _ = recursive_search(chebi_antimetabolite, chebi_has_role, chebi_names, anti, True)
    data, _ = recursive_search(anti, chebi_is_a, chebi_names, data, True)

    data['compound_id'] = data.compound_id.apply(int)
    return data


def parse_pubchem(summary_file):
    pubchem = DataFrame(columns=["name", "molecular_weight", "formula", "uipac_name", "create_date", "compound_id"])

    with open(summary_file) as pubchem_data:
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
    return pubchem
