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

"""
This scripts takes ChEBI relations and vertices to recursive query for entities:
    1. with a "is_a" relationship to any compound marked as analog[ue].
    2. with a "has_role" relationship to any compound marked as antimetabolite (CHEBI:35221).
    2. with a "is_a" relationship to any compound identified in 3.
"""

import os

from marsi.utils import data_dir
from pandas import DataFrame

chebi_names = DataFrame.from_csv(os.path.join(data_dir, "chebi_names_3star.txt"), sep="\t")
chebi_names.fillna("", inplace=True)
chebi_names.index.name = "id"

chebi_names.columns = map(str.lower, chebi_names.columns)
chebi_names.drop_duplicates('compound_id', keep='last', inplace=True)
chebi_names['adapted'] = chebi_names.adapted.apply(lambda v: v == "T")

chebi_analogues = chebi_names[chebi_names.name.str.contains('analog')]
chebi_antimetabolite = chebi_names[chebi_names.compound_id == 35221]

chebi_relations = DataFrame.from_csv(os.path.join(data_dir, "chebi_relation_3star.tsv"), sep="\t")
chebi_relations.columns = map(str.lower, chebi_relations.columns)
chebi_relations.index.name = "id"

chebi_vertices = DataFrame.from_csv(os.path.join(data_dir, "chebi_vertice_3star.tsv"), sep="\t")
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
data.to_csv(os.path.join(data_dir, "chebi_analogues_filtered.csv"))
