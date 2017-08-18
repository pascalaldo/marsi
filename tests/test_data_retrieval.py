# Copyright 2017 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

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
from pandas import DataFrame

from marsi.io import retriaval, parsers

TRAVIS = os.getenv("TRAVIS", False)
if TRAVIS:  # TRAVIS value is 'true'
    TRAVIS = True


@pytest.mark.skipif(TRAVIS, reason="Do not download on travis")
def test_retrieve_bigg(tmpdir):
    bigg_dir = tmpdir.mkdir("bigg")
    dest = bigg_dir.join("bigg_models_reactions.txt")
    retriaval.retrieve_bigg_reactions(dest.strpath)
    statinfo = os.stat(dest.strpath)
    assert statinfo.st_size > 0

    dest = bigg_dir.join("bigg_models_metabolites.txt")
    retriaval.retrieve_bigg_metabolites(dest.strpath)
    statinfo = os.stat(dest.strpath)
    assert statinfo.st_size > 0


@pytest.mark.skipif(TRAVIS, reason="Do not download on travis")
def test_retrieve_drugbank(tmpdir):
    drugbank_dir = tmpdir.mkdir("drugbank")
    dest = drugbank_dir.join("drugbank_open_structures.sdf")
    retriaval.retrieve_drugbank_open_structures(dest=dest.strpath)
    statinfo = os.stat(dest.strpath)
    assert statinfo.st_size > 0

    dest = drugbank_dir.join("drugbank_open_vocabulary.txt")
    retriaval.retrieve_drugbank_open_vocabulary(dest=dest.strpath)
    statinfo = os.stat(dest.strpath)
    assert statinfo.st_size > 0


@pytest.mark.skipif(TRAVIS, reason="Do not download on travis")
def test_retrieve_chebi(tmpdir):
    chebi_dir = tmpdir.mkdir("chebi")
    sdf_dest = chebi_dir.join("chebi_lite_3star.sdf")
    retriaval.retrieve_chebi_structures(dest=sdf_dest.strpath)
    statinfo = os.stat(sdf_dest.strpath)
    assert statinfo.st_size > 0

    names_dest = chebi_dir.join("chebi_names_3star.txt")
    retriaval.retrieve_chebi_names(dest=names_dest.strpath)
    statinfo = os.stat(names_dest.strpath)
    assert statinfo.st_size > 0

    relation_dest = chebi_dir.join("chebi_relation_3star.sdf")
    retriaval.retrieve_chebi_relation(dest=relation_dest.strpath)
    statinfo = os.stat(relation_dest.strpath)
    assert statinfo.st_size > 0

    vertice_dest = chebi_dir.join("chebi_vertice_3star.sdf")
    retriaval.retrieve_chebi_vertice(dest=vertice_dest.strpath)
    statinfo = os.stat(vertice_dest.strpath)
    assert statinfo.st_size > 0

    chebi_data = parsers.parse_chebi_data(names_dest.strpath, vertice_dest.strpath, relation_dest.strpath)
    assert isinstance(chebi_data, DataFrame)
    assert len(chebi_data) > 0
    assert 'compound_id' in chebi_data.columns


@pytest.mark.skipif(TRAVIS, reason="Do not download on travis")
def test_retrieve_brite(tmpdir):
    kegg_dir = tmpdir.mkdir("kegg")
    dest = kegg_dir.join("kegg_brite_08310.keg")
    retriaval.retrieve_kegg_brite(dest=dest.strpath)

    statinfo = os.stat(dest.strpath)
    assert statinfo.st_size > 0

    brite_data = parsers.parse_kegg_brite(dest.strpath)
    assert isinstance(brite_data, DataFrame)

    assert len(brite_data) > 1


@pytest.mark.skip(reason="Takes too long to download, because it is a large file")
def test_retrieve_zinc(tmpdir):
    zinc_dir = tmpdir.mkdir("zinc")
    dest = zinc_dir.join("zinc_16_prop.tsv")
    retriaval.retrieve_zinc_properties(dest=dest.strpath)
    statinfo = os.stat(dest.strpath)
    assert statinfo.st_size > 0

    dest = zinc_dir.join("zinc_16.sdf.gz")
    retriaval.retrieve_zinc_structures(dest=dest.strpath)
    statinfo = os.stat(dest.strpath)
    assert statinfo.st_size > 0
