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
import openbabel
import pybel
import re
from marsi.parallel import apply
from mongoengine import Document, StringField, ListField, ReferenceField, MapField, IntField
from pandas import DataFrame

from marsi.io.plots import summary_plot
from marsi.processing.chemistry import inchi_to_molecule, fingerprint, mol_to_inchi_key, mol_to_inchi

__all__ = ['Database', 'Metabolite']

INCHI_KEY_REGEX = re.compile("[0-9A-Z]{14}\-[0-9A-Z]{8,10}\-[0-9A-Z]")


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

    def __iter__(self):
        return iter(self.collection.objects)

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
    meta = {'collection': 'metabolites'}

    inchi_key = StringField(primary_key=True, regex=INCHI_KEY_REGEX, max_length=27)
    inchi = StringField()
    references = ListField(ReferenceField('Reference'))
    synonyms = ListField(StringField())
    _fingerprints = MapField(ListField(IntField()))

    @classmethod
    def summary(cls, plot=True):
        data = DataFrame(columns=["chebi", "kegg", "pubchem", "drugbank"])
        for metabolite in cls.objects:
            chebi_refs = [r.accession for r in metabolite.references if r.database == "chebi"]
            kegg_refs = [r.accession for r in metabolite.references if r.database == "kegg"]
            pubchem_refs = [r.accession for r in metabolite.references if r.database == "pubchem"]
            drugbank_refs = [r.accession for r in metabolite.references if r.database == "drugbank"]

            if len(chebi_refs) == 0:
                chebi_refs = None
            if len(kegg_refs) == 0:
                kegg_refs = None
            if len(pubchem_refs) == 0:
                pubchem_refs = None
            if len(drugbank_refs) == 0:
                drugbank_refs = None

            data.loc[metabolite.inchi_key] = [chebi_refs, kegg_refs, pubchem_refs, drugbank_refs]

        if plot:
            summary_plot(data)

        return data

    @classmethod
    def apply(cls, func, args=None, context=None, **kwargs):
        if args is None:
            args = tuple()
        if context:
            return [func(getattr(o, context), *args, **kwargs) for o in cls.objects]
        else:
            return apply(func, cls.objects,  *args, **kwargs)

    @classmethod
    def get(cls, inchi_key):
        """
        Retrieves a metabolite using the InChI Key.

        Arguments
        ---------
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

        Arguments
        ---------
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
    def from_molecule(cls, molecule, references):
        try:
            metabolite = cls.get(inchi_key=mol_to_inchi_key(molecule))
            for reference in references:
                if reference not in metabolite.references:
                    metabolite.references.append(reference)
        except KeyError:
            metabolite = cls(inchi_key=mol_to_inchi_key(molecule),
                             inchi=mol_to_inchi(molecule),
                             references=references)
        metabolite.save()

        return metabolite

    def __repr__(self):
        return "Metabolite %s [%s]" % (self.inchi_key, ", ".join(str(r) for r in self.references))

    def __str__(self):
        """
        Metabolites are represented by it's InChI.

        """
        return self.inchi

    def fingerprint(self, fpformat='maccs'):
        if fpformat in self._fingerprints:
            intlist = self._fingerprints[fpformat]
            vector = openbabel.vectorUnsignedInt(len(intlist))
            [vector.__setitem__(i, intlist[i]) for i in range(len(intlist))]
            return pybel.Fingerprint(vector)
        else:
            fp = fingerprint(inchi_to_molecule(str(self)), fpformat=fpformat)
            self._fingerprints[fpformat] = list(fp.fp)
            return fp


class Database(object):
    metabolites = CollectionWrapper(Metabolite)
