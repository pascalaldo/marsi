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

from d3flux import flux_map
from escher import Builder


CENTRAL_CARBON_METABOLISM = ["Citric Acid Cycle", "Glycolysis/Gluconeogenesis"]
DEFAULT_COMPARTMENTS = ['c']


class ModelVisualization(object):
    def __init__(self, model):
        self._model = model
        self._pathways = list(set(reaction.subsystem for reaction in model.reactions))

    def view(self, pathways=None, compartments=DEFAULT_COMPARTMENTS):
        if pathways is None:
            pathways = CENTRAL_CARBON_METABOLISM

        for metabolite in self._model.metabolites:
            metabolite.notes['map_info'] = metabolite.notes.get('map_info', {})
            metabolite.notes['map_info']['hidden'] = any(r.subsystem in pathways for r in metabolite.reactions)

        flux_map(self._model,
                 excluded_compartments=[c for c in self._model.compartments.keys() if c not in compartments])

    @property
    def pathways(self):
        return self._pathways


def view_metabolite_counts(counts, map_name, absolute=True):
    if os.path.exists(map_name):
        map_json = map_name
        map_name = None
    else:
        map_json = None

    if not absolute:
        total = sum(counts.values())
        for m, c in counts:
            counts[m] = c/total

    return Builder(map_name=map_name, map_json=map_json, metabolite_data=counts)
