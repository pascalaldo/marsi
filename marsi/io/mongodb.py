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

import pybel
from cameo.parallel import SequentialView
from mongoengine import Document, StringField, ListField, ReferenceField, MapField, IntField, FloatField, \
    BooleanField, ValidationError
from mongoengine.base import BaseField
from numpy import ndarray, zeros, array
from pandas import DataFrame

from marsi.chemistry import INCHI_KEY_REGEX, openbabel, rdkit
from marsi.io.enrichment import map_uniprot_from_pdb_ids
from marsi.io.plots import summary_plot
from marsi.utils import INCHI_KEY_TYPE, unique

__all__ = ['Database', 'Metabolite']

ENCODING = 'utf-8'


class PointField(BaseField):
    def __init__(self, **kwargs):
        super(PointField, self).__init__(**kwargs)

    # value is (x, y, z)
    def to_python(self, value):
        try:
            value = (value)
        except ValueError:
            pass
        return value

    def validate(self, value, clean=True):
        if len(value) != 3:
            raise ValidationError
        if any(not isinstance(c, float) for c in value):
            raise ValidationError

    def prepare_query_value(self, op, value):
        if value is None:
            return value

        return super(PointField, self).prepare_query_value(op, value)


class ColumnVector(object):
    def __init__(self, collection, column):
        assert issubclass(collection, Document)
        self.collection = collection
        self.column = column

    def apply(self, *args, **kwargs):
        return self.collection.apply(*args, context=self.column, **kwargs)


class CollectionWrapper(object):
    def __init__(self, collection):
        assert issubclass(collection, Document)
        self.collection = collection

    def __len__(self):
        return self.collection.objects.count()

    def __iter__(self):
        return self.collection.objects

    def __getitem__(self, item):
        return self.collection.objects.__getitem__(item)

    def __getattribute__(self, item):
        try:
            return super(CollectionWrapper, self).__getattribute__(item)
        except AttributeError:
            if item in self.collection._fields.keys():
                return ColumnVector(self.collection, item)
            else:
                return getattr(self.collection, item)


class Reference(Document):
    meta = {'collection': 'references'}
    database = StringField(unique_with="accession")
    accession = StringField(unique_with="database")

    @classmethod
    def add_reference(cls, database, accession):
        reference = cls.objects(database=database, accession=accession)
        if reference.count() == 0:
            reference = cls(database=database, accession=accession)
            reference.save()
            return reference
        else:
            return reference.first()

    def __str__(self):
        return "%s: %s" % (self.database, self.accession)


class Metabolite(Document):
    meta = {
        'collection': 'metabolites',
        'index_options': {
            'unique': True,
            'background': True
        },
        'indexes': [
            "#inchi_key"
        ]

    }

    inchi_key = StringField(primary_key=True, regex=INCHI_KEY_REGEX, max_length=27)
    inchi = StringField()
    references = ListField(ReferenceField('Reference'))
    synonyms = ListField(StringField())
    _fingerprints = MapField(ListField(IntField()))
    analog = BooleanField(default=False)
    force_field_coordinates = ListField(PointField())
    sdf = StringField()

    _atoms = IntField(min_value=2)
    _bonds = IntField(min_value=1)

    _volume = FloatField()

    @classmethod
    def summary(cls, plot=True, iterator=iter):
        columns = ["chebi", "kegg", "pubchem", "drugbank", "zinc"]
        rows = cls.objects.count()
        matrix = zeros((rows, len(columns)), dtype=bool)
        index = ndarray((rows, 1), dtype=INCHI_KEY_TYPE)
        for i, metabolite in enumerate(iterator(cls.objects)):
            chebi = False
            kegg = False
            pubchem = False
            drugbank = False
            zinc = False

            for r in metabolite.references:
                chebi = chebi or r.accession == "chebi"
                kegg = kegg or r.database == "kegg"
                pubchem = pubchem or r.database == "pubchem"
                drugbank = drugbank or r.database == "drugbank"
                zinc = zinc or r.database == "zinc15"

            matrix[i, :] = [chebi, kegg, pubchem, drugbank, zinc]
            index[i] = metabolite.inchi_key

        matrix = DataFrame(matrix, index=index[:, 0], columns=columns)

        if plot:
            summary_plot(matrix)

        return matrix

    @classmethod
    def apply(cls, func, args=None, context=None, view=SequentialView(), **kwargs):
        if args is None:
            args = tuple()
        if context:
            return [func(getattr(o, context), *args, **kwargs) for o in cls.objects]
        else:
            return view.apply(func, cls.objects,  *args, **kwargs)

    @classmethod
    def get(cls, inchi_key):
        """
        Retrieves a metabolite using the InChI Key.

        Parameters
        ----------
        inchi_key: str
            A valid InChi Key.

        Returns
        -------
        Metabolite

        Raises
        ------
        KeyError
            If the InChI Key is not available.
        """
        result = cls.objects(inchi_key=inchi_key)
        if result.count() == 0:
            raise KeyError(inchi_key)
        else:
            return result.first()

    @classmethod
    def query(cls, reference):
        """
        Search the metabolites database using references from external database.

        Parameters
        ----------
        reference: str
            A reference.

        Returns
        -------
        list
            A list of Metabolite with the given reference.

        """
        result = cls.objects(references=reference)
        if result.count() == 0:
            raise KeyError(reference)
        else:
            return result.all()

    @classmethod
    def from_molecule(cls, molecule, references, synonyms, analog=False):
        try:
            metabolite = cls.get(inchi_key=openbabel.mol_to_inchi_key(molecule))
            for reference in references:
                if reference not in metabolite.references:
                    metabolite.references.append(reference)
            for synonym in synonyms:
                metabolite.synonyms.append(synonym)
        except KeyError:
            metabolite = cls(inchi_key=openbabel.mol_to_inchi_key(molecule),
                             inchi=openbabel.mol_to_inchi(molecule),
                             references=references,
                             analog=analog,
                             sdf=openbabel.molecule_to_sdf(molecule))
        metabolite.save()

        return metabolite

    @property
    def num_atoms(self):
        if self._atoms is None:
            molecule = self.molecule('openbabel', get3d=False)
            self._atoms = molecule.OBMol.NumAtoms()
            self._bonds = molecule.OBMol.NumBonds()
            self.save()
        return self._atoms

    @property
    def num_bonds(self):
        if self._bonds is None:
            molecule = self.molecule('openbabel', get3d=False)
            self._atoms = molecule.OBMol.NumAtoms()
            self._bonds = molecule.OBMol.NumBonds()
            self.save()
        return self._bonds

    @property
    def atom_coordinates(self):
        if len(self.force_field_coordinates) == 0:
            mol = self.molecule()
            mol.make3D(forcefield='mmff94', steps=int(10e6))
            self.force_field_coordinates = [a.coords for a in mol.atoms]
            self.save()
        return self.force_field_coordinates

    @property
    def formula(self):
        return self.molecule('openbabel', get3d=False).formula


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
            if get3d and self.sdf is not None:
                return rdkit.sdf_to_molecule(self._sdf, from_file=False)
            else:
                return rdkit.inchi_to_molecule(str(self))
        else:
            raise ValueError("Invalid library: %s, please choose between `openbabel` or `rdkit`")

    # NOTE: Hack to get SDF files correct
    @property
    def _sdf(self):
        if self.sdf is None:
            raise ValueError("SDF is not available")
        if self.sdf.startswith("OpenBabel"):
            return "QuickFix1234\n" + self.sdf
        else:
            return self.sdf

    def fingerprint(self, fpformat='ecfp10', hash=False, bits=1024):
        if fpformat not in self._fingerprints:
            fp = openbabel.fingerprint(self.molecule(library='openbabel'), fpformat=fpformat)
            self._fingerprints[fpformat] = list(fp.fp)
            self.save()
        if hash:
            fp = self._fingerprints[fpformat]
            uint_vector = pybel.ob.vectorUnsignedInt(len(fp))
            for i in range(len(fp)):
                uint_vector[i] = fp[i]
            return openbabel.fingerprint_to_bits(pybel.Fingerprint(uint_vector), bits=bits)
        else:
            return array(self._fingerprints[fpformat])

    @property
    def solubility(self):
        molecule = self.molecule(library='openbabel')
        return openbabel.solubility(molecule, log_value=False) * molecule.molwt

    def _repr_html_(self):
        mol = self.molecule(library='openbabel')
        mol.make3D(forcefield='mmff94')
        structure = mol._repr_html_() or openbabel.mol_to_svg(mol)
        return """<table>
    <tbody>
        <tr><td><strong>InChi</strong></td><td>%s</td></tr>
        <tr><td><strong>InChi Key</strong></td><td>%s</td></tr>
        <tr><td><strong>Structure</strong></td><td>%s</td></tr>
        <tr><td><strong>DBs</strong></td><td>%s</td></tr>
        <tr><td><strong>Formula</strong></td><td>%s</td></tr>
    </tbody>
</table>""" % (self.inchi, self.inchi_key, structure, "; ".join(str(r) for r in self.references), mol.formula)

    def __repr__(self):
        return "Metabolite %s [%s]" % (self.inchi_key, ", ".join(str(r) for r in self.references))

    def __str__(self):
        """
        Metabolites are represented by it's InChI.

        """
        return self.inchi

    @property
    def ic50s(self):
        return IC50.objects(metabolite=self)

    @property
    def ec50s(self):
        return EC50.objects(metabolite=self)


class InhibitionMeasurement(Document):
    meta = {
        'collection': 'inhibition_measurement',
        'allow_inheritance': True,
        'index_options': {
            'background': True
        },
        'indexes': [
            'metabolite'
        ]

    }
    value = FloatField(min_value=0)
    operator = StringField(choices=["eq", "gt", "gte", "lt", "lte"])
    metabolite = ReferenceField(Metabolite)
    target_organism = StringField()
    target_name = StringField()
    uniprot_ids = ListField(StringField())

    @property
    def operator_sign(self):
        if self.operator == "eq":
            return "="
        elif self.operator == "gt":
            return ">"
        elif self.operator == "gte":
            return ">="
        elif self.operator == "lt":
            return "<"
        elif self.operator == "lte":
            return "<="

    @classmethod
    def _from_binding_db(cls, row, field, first_time=False):
        metabolite = Metabolite.get(row.inchi_key)

        if isinstance(row[field], float):
            operator, value = "eq", row[field]
        else:
            value = row[field].strip()
            if value.startswith("<"):
                operator, value = "lt", float(value[1:])
            elif value.startswith(">"):
                operator, value = "gt", float(value[1:])
            else:
                operator, value = "eq", float(value)

        entry = None
        if not first_time:
            query = cls.objects(target_name=row.target_name, target_organism=row.target_organism,
                                metabolite=metabolite, value=value, operator=operator)
            if query.count() == 1:
                entry = query[0]

        if first_time or entry is None:
            entry = cls(target_name=row.target_name, target_organism=row.target_organism,
                        metabolite=metabolite, value=value, operator=operator)

        if isinstance(row.pdb_ids, str):
            pdb_ids = tuple(sorted(row.pdb_ids.split(",")))
            uniprot_ids = list(map_uniprot_from_pdb_ids(pdb_ids).values())
            uniprot_ids = sum(uniprot_ids, [])
            uniprot_ids += list(entry.uniprot_ids)
            uniprot_ids = unique(uniprot_ids)
            entry.uniprot_ids = uniprot_ids

        entry.save()
        return entry


class IC50(InhibitionMeasurement):
    """
    Half maximal inhibitory concentration.
    """
    @classmethod
    def ic50_from_binding_db(cls, row, first_time=False):
        return cls._from_binding_db(row, 'ic50', first_time=first_time)

    def __repr__(self):
        return "IC50 %s %f @ (%s %s)" % (self.operator_sign, self.value, self.target_organism, self.target_name)


class EC50(InhibitionMeasurement):
    """
    Half maximal effective concentration.
    """
    @classmethod
    def ec50_from_binding_db(cls, row, first_time=False):
        return cls._from_binding_db(row, 'ec50', first_time=first_time)

    def __repr__(self):
        return "EC50 %s %f @ (%s %s)" % (self.operator_sign, self.value, self.target_organism, self.target_name)


class Database(object):
    metabolites = CollectionWrapper(Metabolite)
    ic50s = CollectionWrapper(IC50)
    ec50s = CollectionWrapper(EC50)
