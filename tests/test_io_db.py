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
import pytest

from openbabel import pybel
import rdkit
from sqlalchemy.exc import IntegrityError

from marsi.chemistry import openbabel

from marsi.io.db import Metabolite, Reference, Database
from marsi.config import default_session


def test_get_metabolite_by_inchi(benchmark):
    met = benchmark(Metabolite.get, "MKUXAQIIEYXACX-UHFFFAOYSA-N", session=default_session)
    assert str(met) == met.inchi

    ob_molecule = met.molecule('openbabel')
    assert isinstance(ob_molecule, pybel.Molecule)

    rd_molecule = met.molecule('rdkit')
    assert isinstance(rd_molecule, rdkit.Chem.rdchem.Mol)

    with pytest.raises(KeyError):
        Metabolite.get("bla-bla-bla")


@pytest.fixture(params=range(5))
def metabolite(request):
    return default_session.query(Metabolite).filter(Metabolite.id == int(request.param) + 1).one()


def _calc_fp(metabolite):
    mol = metabolite.molecule(get3d=False)
    fp = openbabel.fingerprint(mol, 'maccs')
    return openbabel.fingerprint_to_bits(fp, openbabel.fp_bits['maccs'])


def test_fingerprint_method(metabolite, benchmark):
    fp = _calc_fp(metabolite)
    ob_fp = benchmark(_calc_fp, metabolite)
    assert (fp == ob_fp)


def test_collection_wrapper():
    for i in range(10):
        assert Database.metabolites[i] == default_session.query(Metabolite).filter(Metabolite.id == i + 1).one()


def test_add_reference():
    database = "test_db"
    accession1 = "entry1"
    ref1 = Reference.add_reference(database, accession1, default_session)

    assert ref1.database == database
    assert ref1.accession == accession1

    accession2 = "entry2"

    ref2 = Reference.add_reference(database, accession2)

    assert ref2.database == database
    assert ref2.accession == accession2

    ref3 = Reference.add_reference(database, accession2)

    assert ref3.database == database
    assert ref3.accession == accession2

    default_session.commit()

    with pytest.raises(IntegrityError):
        ref4 = Reference(database=database, accession=accession1)
        default_session.add(ref4)
        default_session.commit()

    default_session.rollback()

    default_session.delete(ref1)
    default_session.delete(ref2)
    default_session.delete(ref3)
    default_session.commit()
