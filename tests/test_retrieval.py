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
import unittest

import os
import tempfile

from marsi.io.retriaval import retrieve_bigg_metabolites, retrieve_bigg_reactions
from marsi.io.retriaval import retrieve_chebi_names, retrieve_chebi_relation, retrieve_chebi_structures,\
    retrieve_chebi_vertice
from marsi.io.retriaval import retrieve_drugbank_open_structures, retrieve_drugbank_open_vocabulary
from marsi.io.retriaval import retrieve_kegg_brite


class RetrievalTestCase(unittest.TestCase):
    temp_dir = None

    @classmethod
    def setUpClass(cls):
        cls.temp_dir = tempfile.gettempdir()

    def assertFileExistsAndHasContent(self, filepath):
        self.assertTrue(os.path.isfile(filepath))
        self.assertTrue(os.path.exists(filepath))
        self._head(filepath)
        self.assertGreater(os.path.getsize(filepath), 0)

    @staticmethod
    def _head(filepath):
        with open(filepath) as file:
            print("FILE: %s" % filepath)
            for i in range(10):
                line = next(file)
                if len(line) > 150:
                    print("%i: %s" % (i+1, line[:150] + " ..."))
                else:
                    print("%i: %s" % (i+1, line))

    def test_retrieve_bigg(self):
        metabolites = os.path.join(self.temp_dir, "bigg_metabolites.txt")
        retrieve_bigg_metabolites(metabolites)
        self.assertFileExistsAndHasContent(metabolites)
        os.remove(metabolites)

        reactions = os.path.join(self.temp_dir, "bigg_reactions.txt")
        retrieve_bigg_reactions(reactions)
        self.assertFileExistsAndHasContent(reactions)
        os.remove(reactions)

    def test_retrieve_chebi(self):
        names = os.path.join(self.temp_dir, "chebi_names")
        retrieve_chebi_names(names)
        self.assertFileExistsAndHasContent(names)
        os.remove(names)

        relation = os.path.join(self.temp_dir, "chebi_relation")
        retrieve_chebi_relation(relation)
        self.assertFileExistsAndHasContent(relation)
        os.remove(relation)

        structures = os.path.join(self.temp_dir, "chebi_structures")
        retrieve_chebi_structures(structures)
        self.assertFileExistsAndHasContent(structures)
        os.remove(structures)

        vertice = os.path.join(self.temp_dir, "chebi_vertice")
        retrieve_chebi_vertice(vertice)
        self.assertFileExistsAndHasContent(vertice)
        os.remove(vertice)

    def test_retrieve_kegg(self):
        kegg_brite = os.path.join(self.temp_dir, "kegg_brite")
        retrieve_kegg_brite(kegg_brite)
        self.assertFileExistsAndHasContent(kegg_brite)
        os.remove(kegg_brite)

    def test_retrieve_drugbank(self):
        structures = os.path.join(self.temp_dir, "drugbank_strucutres")
        retrieve_drugbank_open_structures(dest=structures)
        self.assertFileExistsAndHasContent(structures)
        os.remove(structures)

        vocabulary = os.path.join(self.temp_dir, "drugbank_vocabulary")
        retrieve_drugbank_open_vocabulary(dest=vocabulary)
        self.assertFileExistsAndHasContent(vocabulary)
        os.remove(vocabulary)


