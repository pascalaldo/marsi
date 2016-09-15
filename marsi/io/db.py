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


import re

from mongoengine import Document, StringField, ListField
from marsi.parallel import apply

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

    def __getattribute__(self, item):
        try:
            return super(CollectionWrapper, self).__getattribute__(item)
        except AttributeError:
            if item in self.collection._fields.keys():
                return ColumnVector(self.collection, item)
            else:
                raise AttributeError(item)


class Metabolite(Document):

    meta = {'collection': 'metabolites'}

    inchi_key = StringField(primary_key=True, regex=INCHI_KEY_REGEX, max_length=27)
    inchi = StringField()
    references = ListField(StringField())

    @classmethod
    def apply(cls, func, args=None, context=None, **kwargs):
        if args is None:
            args = tuple()
        if context:
            return apply([func(getattr(o, context), *args, **kwargs) for o in cls.objects])
        else:
            return apply([func(o, *args, **kwargs) for o in cls.objects])

    @classmethod
    def get(cls, inchi_key):
        result = cls.objects(inchi_key=inchi_key)
        if result.count() == 0:
            raise KeyError(inchi_key)
        else:
            return result.first()

    @classmethod
    def query(cls, reference):
        result = cls.objects(references=reference)
        if result.count() == 0:
            raise KeyError(reference)
        else:
            return result.all

    def __repr__(self):
        return "Metabolite %s [%s]" % (self.inchi_key, ", ".join(self.references))


class Database(object):
    metabolites = CollectionWrapper(Metabolite)