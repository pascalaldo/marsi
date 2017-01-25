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
from functools import lru_cache
from pandas import DataFrame

from marsi import bigg_api
from marsi.cobra.flux_analysis import compete_metabolite, inhibit_metabolite
from marsi.io.bigg import bigg_metabolites
from marsi.algorithms.enrichment import inchi_from_chebi, inchi_from_kegg

__all__ = ['find_inchi_for_bigg_metabolite', 'degree', 'in_degree', 'out_degree']

logger = logging.getLogger(__name__)

DATABASE_LINKS = 'database_links'

CHEBI = 'CHEBI'
KEGG = "KEGG Compound"


@lru_cache(maxsize=1024)
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


def degree(metabolite):
    """
    Degree of a metabolite is the number of reactions coming in and out.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A Metabolite.

    Returns
    -------
    int
        The degree.

    """
    return sum(2 if r.reversibility else 1 for r in metabolite.reactions)


def in_degree(metabolite):
    """
    In degree of a metabolite is the number of reactions coming in.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A Metabolite.

    Returns
    -------
    int
        The degree.

    """
    return len([r for r in metabolite.reactions if r.metabolites[metabolite] < 0 or r.reversibility])


def out_degree(metabolite):
    """
    Degree of a metabolite is the number of reactions coming out.

    Parameters
    ----------
    metabolite: cobra.Metabolite
        A Metabolite.

    Returns
    -------
    int
        The degree.

    """
    return len([r for r in metabolite.reactions if r.metabolites[metabolite] > 0 or r.reversibility])


def apply_anti_metabolite(metabolites, essential_metabolites, reference, inhibition_fraction=.0,
                          competition_fraction=.0, allow_accumulation=True, ignore_transport=True, time_machine=None):
    exchanges = set()

    if any(met.id in essential_metabolites for met in metabolites):
        for metabolite in metabolites:
            exchanges.add(compete_metabolite(metabolite,
                                             reference,
                                             allow_accumulation=allow_accumulation,
                                             ignore_transport=ignore_transport,
                                             fraction=competition_fraction,
                                             time_machine=time_machine))
    else:
        for metabolite in metabolites:
            exchanges.add(inhibit_metabolite(metabolite,
                                             reference,
                                             allow_accumulation=allow_accumulation,
                                             ignore_transport=ignore_transport,
                                             fraction=inhibition_fraction,
                                             time_machine=time_machine))

    return exchanges





