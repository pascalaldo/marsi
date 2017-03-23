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

import pytest
import requests
import six

from marsi import bigg_api

METABOLITES = {
    "glc__D": dict(name="D-Glucose", bigg_id="glc__D"),
    "ser__L": dict(name="L-Serine", bigg_id="ser__L"),
    "g3p": dict(name="Glyceraldehyde 3-phosphate", bigg_id="g3p")
}


REACTIONS = {
    "GAPD": dict(name="Glyceraldehyde-3-phosphate dehydrogenase", bigg_id="GAPD", pseudoreaction=False),
    "ADA": dict(name="Adenosine deaminase", bigg_id="ADA", pseudoreaction=False)
}


@pytest.fixture(params=list(METABOLITES.keys()))
def metabolite(request):
    return METABOLITES[request.param]


@pytest.fixture(params=list(REACTIONS.keys()))
def reaction(request):
    return REACTIONS[request.param]


@pytest.fixture(params=["xml", "xml-gz", "mat", "json"])
def model_format(request):
    return request.param.replace("-", ".")


@pytest.fixture
def bigg_model():
    return {
        "bigg_id": "iND750",
        "organism": "Saccharomyces cerevisiae S288c",
        "metabolite_count": 1059,
        "gene_count": 750,
        "reaction_count": 1266,
        "reference_type": "pmid",
        "reference_id": "15197165"
    }


def test_get_database_version():
    db_version = bigg_api.database_version()
    assert isinstance(db_version, bigg_api.DBVersion)
    assert db_version.api_version == "v2"

    assert str(db_version).startswith("BiGG Database version")


def test_model_details(bigg_model):
    yeast_details = bigg_api.model_details(bigg_model['bigg_id'])
    assert yeast_details["organism"] == bigg_model["organism"]
    assert yeast_details["metabolite_count"] == bigg_model["metabolite_count"]
    assert yeast_details["gene_count"] == bigg_model["gene_count"]
    assert yeast_details["reaction_count"] == bigg_model["reaction_count"]
    assert yeast_details["reference_type"] == bigg_model["reference_type"]
    assert yeast_details["reference_id"] == bigg_model["reference_id"]


def test_download_model(model_format, bigg_model, tmpdir):
    path = tmpdir.mkdir("models")
    final_doc = path.join("%s.%s" % (bigg_model['bigg_id'], model_format))
    bigg_api.download_model(bigg_model['bigg_id'], model_format, True, path.strpath)
    assert os.path.exists(final_doc.strpath)
    assert os.path.isfile(final_doc.strpath)


def test_model_reactions(bigg_model):
    result = bigg_api.list_model_reactions(bigg_model['bigg_id'])
    assert result['results_count'] == bigg_model['reaction_count']

    for i in range(10):
        reaction_info = result['results'][i]
        reaction_result = bigg_api.get_model_reaction(bigg_model['bigg_id'], reaction_info['bigg_id'])
        assert reaction_result['bigg_id'] == reaction_info['bigg_id']
        assert reaction_result['name'] == reaction_info['name']

    with pytest.raises(requests.HTTPError):
        bigg_api.get_model_reaction(bigg_model['bigg_id'], "non-existent")


def test_model_metabolites(bigg_model):
    result = bigg_api.list_model_metabolites(bigg_model['bigg_id'])
    assert result['results_count'] == bigg_model['metabolite_count']

    for i in range(10):
        metabolite_info = result['results'][i]
        metabolite_id = metabolite_info['bigg_id'] + "_" + metabolite_info['compartment_bigg_id']
        metabolite_result = bigg_api.get_model_metabolite(bigg_model['bigg_id'], metabolite_id)
        assert metabolite_result['bigg_id'] == metabolite_info['bigg_id']
        assert metabolite_result['name'] == metabolite_info['name']

    with pytest.raises(requests.HTTPError):
        bigg_api.get_model_metabolite(bigg_model['bigg_id'], "non-existent")


def test_get_metabolite(metabolite):
    result = bigg_api.get_metabolite(metabolite['bigg_id'])
    for key, value in six.iteritems(metabolite):
        assert result[key] == value


def test_fail_get_metabolite():
    with pytest.raises(requests.HTTPError):
        bigg_api.get_metabolite("non-existent")


def test_get_reaction(reaction):
    result = bigg_api.get_reaction(reaction['bigg_id'])
    for key, value in six.iteritems(reaction):
        assert result[key] == value


def test_fail_get_reaction():
    with pytest.raises(requests.HTTPError):
        bigg_api.get_reaction("non-existent")


def test_list_reactions():
    reactions = bigg_api.list_reactions()
    assert 'results_count' in reactions
    assert 'results' in reactions
    assert len(reactions['results']) == reactions['results_count']

    for i in range(1000):
        assert 'bigg_id' in reactions['results'][i]
        assert 'name' in reactions['results'][i]


def test_list_metabolites():
    metabolites = bigg_api.list_metabolites()
    assert 'results_count' in metabolites
    assert 'results' in metabolites
    assert len(metabolites['results']) == metabolites['results_count']

    for i in range(1000):
        assert 'bigg_id' in metabolites['results'][i]
        assert 'name' in metabolites['results'][i]
