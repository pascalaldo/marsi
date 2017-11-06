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
import six
from bitarray import bitarray
from sqlalchemy import func
from sqlalchemy import inspect

from marsi.chemistry import rdkit

from marsi.chemistry import openbabel

from sqlalchemy import Boolean, Integer, String, Table, Text
from sqlalchemy import Column, ForeignKey, Index, UniqueConstraint
from sqlalchemy import TypeDecorator
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import validates, relationship, backref
from sqlalchemy.orm.collections import attribute_mapped_collection
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.sql.elements import and_

from marsi.chemistry.common import INCHI_KEY_REGEX
from marsi.config import default_session

__all__ = ['Database', 'Metabolite', 'Reference']


Base = declarative_base()


class ColumnVector(object):
    def __init__(self, collection, session, column):
        assert issubclass(collection, Base)
        self.collection = collection
        self.session = session
        self.column = column

    def apply(self, *args, **kwargs):
        return self.collection.apply(*args, context=self.column, **kwargs)


class CollectionWrapper(object):
    def __init__(self, collection, session=default_session):
        assert issubclass(collection, Base)
        self.collection = collection
        self.session = session
        self.mapper = inspect(self.collection)

    def dump(self, i):
        return [self.__getitem__(item).dump() for item in range(i)]

    def restore(self, dump, session=default_session):
        for _dump in dump:
            self.collection.restore(_dump, session=session)

    def __len__(self):
        return self.session.query(func.count(self.collection.id))

    def __iter__(self):
        return self.session.query(self.collection).yield_per(1000)

    def __getitem__(self, item):
        return self.session.query(self.collection).filter(self.collection.id == item + 1).one()

    def __getattribute__(self, item):
        try:
            return super(CollectionWrapper, self).__getattribute__(item)
        except AttributeError:
            if item in [attr.key for attr in self.mapper.attrs]:
                return ColumnVector(self.collection, self.session, item)
            else:
                return getattr(self.collection, item)


class Fingerprint(TypeDecorator):
    impl = String

    def process_bind_param(self, value, dialect):
        return value.to01()

    def process_result_value(self, value, dialect):
        return bitarray(value)

    def copy(self, **kw):
        return Fingerprint(self.impl.length)


references_table = Table('metabolite_references', Base.metadata,
                         Column('metabolite_id', Integer, ForeignKey('metabolites.id')),
                         Column('reference_id', Integer, ForeignKey('references.id')))

synonyms_table = Table('metabolite_synonyms', Base.metadata,
                       Column('metabolite_id', Integer, ForeignKey('metabolites.id')),
                       Column('synonym_id', Integer, ForeignKey('synonyms.id')))


class Synonym(Base):
    __tablename__ = "synonyms"
    id = Column(Integer, primary_key=True)
    synonym = Column(String(500), nullable=False)

    __table_args__ = (Index('uq_synonyms', synonym, unique=True), )

    @classmethod
    def add_synonym(cls, synonym, session=default_session):
        query = session.query(cls).filter(cls.synonym == synonym)

        try:
            synonym = query.one()
        except NoResultFound:
            synonym = cls(synonym=synonym)
            session.add(synonym)
            session.flush()

        return synonym

    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns}


class Reference(Base):
    __tablename__ = "references"

    id = Column(Integer, primary_key=True)
    database = Column(String(100), nullable=False)
    accession = Column(String(100), nullable=False)

    __table_args__ = (UniqueConstraint('database', 'accession', name='_database_accession_uc'), )

    @classmethod
    def get(cls, database, accession, session=default_session):
        database = database.strip()
        accession = accession.strip()
        query = session.query(cls).filter(and_(cls.database == database, cls.accession == accession))

        return query.one()

    @classmethod
    def add_reference(cls, database, accession, session=default_session):
        database = database.strip()
        accession = accession.strip()
        query = session.query(cls).filter(and_(cls.database == database, cls.accession == accession))

        try:
            reference = query.one()
        except NoResultFound:
            reference = cls(database=database, accession=accession)
            session.add(reference)
            session.flush()

        return reference

    def __str__(self):
        return "%s: %s" % (self.database, self.accession)

    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns}


class MetaboliteFingerprint(Base):
    __tablename__ = 'metabolite_fingerprints'

    id = Column(Integer, primary_key=True)
    metabolite_id = Column(Integer, ForeignKey('metabolites.id'))
    fingerprint_type = Column(String(10), nullable=False)
    fingerprint = Column(Fingerprint(2048), nullable=False)

    metabolite = relationship("Metabolite", backref=backref(
        "_fingerprints",
        collection_class=attribute_mapped_collection("fingerprint_type"),
        cascade="all, delete-orphan"))

    __table_args__ = (
        UniqueConstraint('metabolite_id', 'fingerprint_type', name='_fp_type_uc'),
    )

    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns}


class Metabolite(Base):
    __tablename__ = "metabolites"

    id = Column(Integer, primary_key=True)
    inchi_key = Column(String(27), nullable=False)
    inchi = Column(String(5000), nullable=False)
    analog = Column(Boolean, default=False)
    formula = Column(String(500), nullable=False)
    num_atoms = Column(Integer, nullable=False)
    num_bonds = Column(Integer, nullable=False)
    num_rings = Column(Integer, nullable=False)
    sdf = Column(Text, nullable=True)

    # solubility = Column(Float)

    references = relationship("Reference", secondary=references_table)
    synonyms = relationship("Synonym", secondary=synonyms_table)
    fingerprints = association_proxy('_fingerprints', 'fingerprint',
                                     creator=lambda key, value: MetaboliteFingerprint(fingerprint_type=key,
                                                                                      fingerprint=value))

    __table_args__ = (
        Index('uq_inchi_key', inchi_key, unique=True),
    )

    @validates('inchi_key')
    def validate_inchi_key(self, key, inchi_key):
        if not INCHI_KEY_REGEX.match(inchi_key):
            raise ValueError("InChI Key %s is not valid" % inchi_key)
        else:
            return inchi_key

    @classmethod
    def get(cls, inchi_key, session=default_session):
        """
        Retrieves a metabolite using the InChI Key.

        Parameters
        ----------
        inchi_key : str
            A valid InChi Key.
        session : sqlalchemy.orm.session.Session
            A database session.

        Returns
        -------
        Metabolite

        Raises
        ------
        KeyError
            If the InChI Key is not available.
        """
        query = session.query(cls).filter(cls.inchi_key == inchi_key)
        try:
            return query.one()
        except NoResultFound:
            raise KeyError(inchi_key)

    @classmethod
    def from_references(cls, references, session=default_session):
        hits = []
        for r in references:
            query = session.query(cls).join(cls.references).filter(Reference.id == r.id)
            hits.extend(list(query.all()))

        return hits

    @classmethod
    def from_molecule(cls, molecule, references, synonyms, analog=False, session=default_session, first_time=False):
        inchi_key = openbabel.mol_to_inchi_key(molecule)
        try:
            if first_time:
                raise KeyError()

            metabolite = cls.get(inchi_key=inchi_key)
            for reference in references:
                if reference not in metabolite.references:
                    metabolite.references.append(reference)
            for synonym in synonyms:
                if synonym not in metabolite.synonyms:
                    metabolite.synonyms.append(synonym)
        except KeyError:
            metabolite = Metabolite(inchi_key=inchi_key,
                                    inchi=openbabel.mol_to_inchi(molecule),
                                    analog=analog,
                                    formula=molecule.formula,
                                    sdf=openbabel.molecule_to_sdf(molecule),
                                    num_atoms=molecule.OBMol.NumAtoms(),
                                    num_bonds=molecule.OBMol.NumBonds(),
                                    num_rings=len(molecule.OBMol.GetSSSR()))
            for reference in references:
                if reference not in metabolite.references:
                    metabolite.references.append(reference)
            for synonym in synonyms:
                if synonym not in metabolite.synonyms:
                    metabolite.synonyms.append(synonym)

            fingerprint = openbabel.fingerprint(molecule, 'maccs')
            bits = openbabel.fp_bits.get('maccs', 2048)
            metabolite.fingerprints['maccs'] = openbabel.fingerprint_to_bits(fingerprint, bits)

            session.add(metabolite)

        return metabolite

    # NOTE: Hack to get SDF files correct
    @property
    def _sdf(self):
        if self.sdf is None:
            raise ValueError("SDF is not available")
        if self.sdf.startswith("OpenBabel"):
            return "QuickFix1234\n" + self.sdf
        else:
            return self.sdf

    def calc_solubility(self):
        molecule = self.molecule(library='openbabel')
        return openbabel.solubility(molecule, log_value=False) * molecule.molwt

    @property
    def volume(self):
        mol = self.molecule(library='openbabel')
        return openbabel.monte_carlo_volume(mol, self.atom_coordinates, tolerance=1, max_iterations=100)

    def molecule(self, library='openbabel', get3d=True):
        if library == 'openbabel':
            if get3d and self.sdf is not None:
                molecule = openbabel.sdf_to_molecule(self._sdf, from_file=False)
                molecule.title = ""
                return molecule
            else:
                return openbabel.inchi_to_molecule(str(self))
        elif library == 'rdkit':
            try:
                if get3d and self.sdf is not None:
                    return rdkit.sdf_to_molecule(self._sdf, from_file=False)
                else:
                    return rdkit.inchi_to_molecule(str(self))
            except Exception as e:
                print(self.inchi_key)
                raise e
        else:
            raise ValueError("Invalid library: %s, please choose between `openbabel` or `rdkit`")

    def fingerprint(self, fingerprint_format='maccs'):
        if fingerprint_format not in self.fingerprints:
            ob_molecule = self.molecule(get3d=False)
            fingerprint = openbabel.fingerprint(ob_molecule, fingerprint_format)
            bits = openbabel.fp_bits.get(fingerprint_format, 2048)
            self.fingerprints[fingerprint_format] = openbabel.fingerprint_to_bits(fingerprint, bits)

        return self.fingerprints[fingerprint_format]

    def _repr_html_(self):
        mol = self.molecule(library='openbabel')
        mol.make3D(forcefield='mmff94')
        mol.removeh()
        structure = mol._repr_html_() or openbabel.mol_to_svg(mol)
        references = "; ".join(str(r) for r in self.references)
        synonyms = ", ".join([s.synonym for s in self.synonyms])
        return """
<table>
    <tbody>
        <tr><td><strong>Formula</strong></td><td>%s</td></tr>
        <tr><td><strong>InChi</strong></td><td>%s</td></tr>
        <tr><td><strong>InChi Key</strong></td><td>%s</td></tr>
        <tr><td><strong>Structure</strong></td><td>%s</td></tr>
        <tr><td><strong>DBs</strong></td><td>%s</td></tr>
        <tr><td><strong>Synonyms</strong></td><td>%s</td></tr>
    </tbody>
</table>""" % (mol.formula, self.inchi, self.inchi_key, structure, references, synonyms)

    def __repr__(self):
        return "Metabolite %s [%s]" % (self.inchi_key, ", ".join(str(r) for r in self.references))

    def __str__(self):
        """
        Metabolites are represented by it's InChI.

        """
        return self.inchi

    def to_dict(self):
        return {c.name: getattr(self, c.name) for c in self.__table__.columns}

    def dump(self):
        references = [ref.to_dict() for ref in self.references]
        synonyms = [syn.to_dict() for syn in self.synonyms]
        fingerprints = {key: fp.to01() for key, fp in six.iteritems(self.fingerprints)}
        metabolite = self.to_dict()

        return dict(metabolite=metabolite, references=references,
                    synonyms=synonyms, fingerprints=fingerprints)

    @classmethod
    def restore(cls, dump, session=default_session):
        metabolite = cls(**dump['metabolite'])
        references = [Reference.add_reference(ref['database'], ref['accession'], session) for ref in dump['references']]
        synonyms = [Synonym.add_synonym(syn['synonym'], session) for syn in dump['synonyms']]
        metabolite.references = references
        metabolite.synonyms = synonyms
        for key, fingerprint in six.iteritems(dump['fingerprints']):
            metabolite.fingerprints[key] = bitarray(fingerprint)
        session.add(metabolite)
        session.commit()


class Database:
    metabolites = CollectionWrapper(Metabolite)
