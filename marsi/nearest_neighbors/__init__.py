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
import math
import os

import numpy as np
from mongoengine import connect

from cameo.parallel import SequentialView

from marsi import config
from marsi.nearest_neighbors.model import NearestNeighbors, DistributedNearestNeighbors

from marsi.utils import data_dir, INCHI_KEY_TYPE, unpickle_large, pickle_large
from marsi.io.mongodb import Database
from marsi.chemistry import SOLUBILITY

__all__ = ['build_nearest_neighbors_model', 'load_nearest_neighbors_model']

MODEL_FILE = os.path.join(data_dir, "fingerprints_default_%s_sol_%s.pickle")


class FeatureReader(object):
    """
    Accessory class to build fingerprints in chunks.

    Attributes
    ----------
    db : str
        The name of the database to connect to.
    fpformat : str
        The fingerprint format (see pybel.fps).
    solubility : str
        One of 'high', 'medium', 'low', 'all'.
    connection_args : dict
        Other arguments of mongoengine.connect


    """
    def __init__(self, db, fpformat='maccs', solubility='high', **connection_args):
        self.db = db
        self.fpformat = fpformat
        self.solubility = solubility
        self.connection_args = connection_args

    def __call__(self, index):
        connect(self.db, **self.connection_args)
        subset = Database.metabolites[index[0]:index[1]]
        indices = []
        fingerprints = []
        fingerprint_lengths = []
        for m in subset:
            if SOLUBILITY[self.solubility](m.solubility):
                fingerprint = m.fingerprint(fpformat=self.fpformat)
                fingerprints.append(fingerprint)
                indices.append(m.inchi_key)
                fingerprint_lengths.append(len(fingerprint))

        _indices = np.ndarray((len(indices), 1), dtype=INCHI_KEY_TYPE)
        for i in range(_indices.shape[0]):
            _indices[i] = indices[i]
        del indices
        return _indices, fingerprints, fingerprint_lengths


def _build_feature_table(database, fpformat='ecfp10', chunk_size=None, solubility='high',
                         database_name=config.db_name, view=SequentialView()):
    reader = FeatureReader(database_name, fpformat=fpformat, solubility=solubility)
    chunk_size = math.ceil(chunk_size)
    n_chunks = math.ceil(len(database)/chunk_size)
    chunks = [((i-1) * chunk_size, i*chunk_size) for i in range(1, n_chunks + 1)]
    res = view.map(reader, chunks)
    indices = np.ndarray((0, 1), dtype=INCHI_KEY_TYPE)
    fingerprints = []
    fingerprint_lengths = []
    for r in res:
        indices = np.concatenate([indices, r[0]])
        fingerprints += r[1]
        fingerprint_lengths += r[2]
    return indices, fingerprints, fingerprint_lengths


def _build_nearest_neighbors_model(indices, features, lengths, n_models):
    chunk_size = math.ceil(len(indices)/n_models)
    chunks = [((i-1) * chunk_size, i*chunk_size) for i in range(1, n_models+1)]
    models = []
    for start, end in chunks:
        models.append(NearestNeighbors(indices[start:end], features[start:end], lengths[start:end]))

    return DistributedNearestNeighbors(models)


def build_nearest_neighbors_model(database, fpformat='fp4', solubility='high', n_models=5,
                                  chunk_size=1e6, view=SequentialView):
    """
    Loads a NN model.

    If a 'default_model.pickle' exists in data it will load the model. Otherwise it will build a model from the
    Database. This can take several hours depending on the size of the database.

    Parameters
    ----------
    database : marsi.io.mongodb.CollectionWrapper
        A Database interface to the metabolites.
    chunk_size : int
        Maximum number of entries per chunk.
    fpformat : str
        The format of the fingerprint (see pybel.fps)
    solubility : str
        One of high, medium, low or all.
    view : cameo.parallel.SequentialView, cameo.parallel.MultiprocesingView
        A view to control parallelization.
    n_models : int
        The number of NearestNeighbors models.
    """

    indices, features, lens = _build_feature_table(database, fpformat=fpformat, chunk_size=chunk_size,
                                                   solubility=solubility, view=view)
    return _build_nearest_neighbors_model(indices, features, lens, n_models)


def load_nearest_neighbors_model(chunk_size=1e6, fpformat="fp4", solubility='all',
                                 view=SequentialView(), model_size=100000):
    """
    Loads a NN model.

    If a 'default_model.pickle' exists in data it will load the model. Otherwise it will build a model from the
    Database. This can take several hours depending on the size of the database.

    Parameters
    ----------
    chunk_size : int
        Maximum number of entries per chunk.
    fpformat : str
        The format of the fingerprint (see pybel.fps)
    solubility : str
        One of high, medium, low or all.
    view : cameo.parallel.SequentialView, cameo.parallel.MultiprocesingView
        A view to control parallelization.
    model_size : int
        The size of each NearestNeighbor in the ensemble.
    """
    if solubility not in SOLUBILITY:
        raise ValueError('%s not one of %s' % (solubility, ", ".join(SOLUBILITY.keys())))

    model_file = MODEL_FILE % (fpformat, solubility)
    if os.path.exists(model_file):
        _indices, _features, _lengths = unpickle_large(model_file, progress=True)
    else:
        print("Building search model (fp: %s, solubility: %s)" % (fpformat, solubility))
        _indices, _features, _lengths = _build_feature_table(Database.metabolites,
                                                             chunk_size=chunk_size,
                                                             fpformat=fpformat,
                                                             solubility=solubility,
                                                             view=view)
        pickle_large((_indices, _features, _lengths), model_file, progress=True)

    n_models = math.ceil(len(_indices) / model_size)
    nn_model = _build_nearest_neighbors_model(_indices, _features, _lengths, n_models)
    return nn_model
