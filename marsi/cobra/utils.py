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
import logging

from collections import Counter

from IProgress.progressbar import ProgressBar
from IProgress.widgets import Percentage, Bar, ETA

from cachetools import cached, LRUCache

from marsi import bigg_api
from marsi.io.bigg import bigg_metabolites
from marsi.io.enrichment import inchi_from_chebi, inchi_from_kegg


__all__ = ['find_inchi_for_bigg_metabolite', 'annotate_metabolite', 'annotate_model']


lru_cache = LRUCache(maxsize=1024)

logger = logging.getLogger(__name__)

CURRENCY_METABOLITES = ["atp", "adp", "nad", "nadh", "nadp", "nadph", "amp",
                        "h2o", "h", "coa", "acp", "pi", 'pppi', 'ppi']

DATABASE_LINKS = 'database_links'

CHEBI = 'CHEBI'
KEGG = "KEGG Compound"


@cached(lru_cache)
def find_inchi_for_bigg_metabolite(model_id, metabolite_id):
    try:
        links = bigg_metabolites.loc[metabolite_id].database_links
    except KeyError:
        metabolite_data = bigg_api.get_model_metabolite(model_id, metabolite_id)
        links = metabolite_data[DATABASE_LINKS]
    inchi_keys = []
    if CHEBI in links:
        inchi_keys += [inchi_from_chebi(link['id']) for link in links[CHEBI]]

    if KEGG in links:
        inchi_keys += [inchi_from_kegg(link['id']) for link in links[KEGG]]

    counter = Counter(inchi_keys)
    del counter[None]

    if len(counter) > 0:
        return max(counter, key=lambda key: counter[key])
    else:
        raise ValueError(metabolite_id)


def annotate_metabolite(metabolite):
    try:
        metabolite.annotation['inchi']
    except KeyError:
        try:
            metabolite.annotation['inchi'] = find_inchi_for_bigg_metabolite(metabolite.model, metabolite.id)
        except ValueError:
            pass


def annotate_model(model):
    pbar = ProgressBar(maxval=len(model.metabolites), widgets=["Annotating: ", Percentage(), Bar(), ETA()])
    [annotate_metabolite(m) for m in pbar(model.metabolites)]


def essential_species_ids(model):
    return {m.id[:-2] for m in model.essential_metabolites()}
