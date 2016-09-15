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
from collections import defaultdict

from marsi.io.db import Database
from marsi.design.operations import NearestNeighbors

from marsi.parallel import apply


class GenericDesignMethod(object):

    def __init__(self, model=None, metabolite_database=Database, min_tanimoto=0.8, min_rmsd=0.8):
        self.model = model
        self.database = metabolite_database
        self.min_tanimoto = min_tanimoto
        self.min_rmsd = min_rmsd
        self.metabolite_neighbors = defaultdict(lambda: [])
        self.map_model()

    def map_model(self):
        nn_function = NearestNeighbors(self.database)
        neighbors = apply(nn_function, self.model.metabolites, min_tanimoto=self.min_tanimoto, min_rmsd=self.min_rmsd)
        self.metabolite_neighbors.update({m.id: n[2] for m, n in zip(self.model.metabolites, neighbors)})

    def __call__(self, target, maximize=True, **kwargs):
        raise NotImplementedError


class RandomMutagenesisDesign(GenericDesignMethod):
    """
    Apply only knockout like designs where total loss of functions are expected.

    """
    def __call__(self, target, maximize=True, **kwargs):
        pass


class ALEDesign(GenericDesignMethod):
    """
    Apply both knockout and flux modulation. The strains will be selected amongst the fastest growers.


    """
    def __call__(self, target, maximize=True, **kwargs):
        pass
