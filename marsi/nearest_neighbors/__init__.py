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
import multiprocessing
from queue import Empty

import numpy as np
import six
from IProgress import ProgressBar, Bar, ETA
from mongoengine import connect

from cameo.parallel import SequentialView
from pandas import DataFrame
from sqlalchemy import create_engine
from sqlalchemy.orm import sessionmaker

from marsi import config
from marsi.chemistry import rdkit
from marsi.chemistry.molecule import Molecule
from marsi.config import default_session, connection_pool
from marsi.io.db import Metabolite
from marsi.nearest_neighbors.model import NearestNeighbors, DistributedNearestNeighbors, DBNearestNeighbors

from marsi.utils import data_dir, INCHI_KEY_TYPE, unpickle_large, pickle_large
from marsi.io.mongodb import Database
from marsi.chemistry import SOLUBILITY

from sqlalchemy import and_


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


def build_feature_table(database, fpformat='ecfp10', chunk_size=None, solubility='high',
                        database_name=config.db_name, view=SequentialView()):
    reader = FeatureReader(database_name, fpformat=fpformat, solubility=solubility)
    chunk_size = math.ceil(chunk_size)
    n_chunks = math.ceil(len(database) / chunk_size)
    chunks = [((i - 1) * chunk_size, i * chunk_size) for i in range(1, n_chunks + 1)]
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
    chunk_size = math.ceil(len(indices) / n_models)
    chunks = [((i - 1) * chunk_size, i * chunk_size) for i in range(1, n_models + 1)]
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

    indices, features, lens = build_feature_table(database, fpformat=fpformat, chunk_size=chunk_size,
                                                  solubility=solubility, view=view)
    return _build_nearest_neighbors_model(indices, features, lens, n_models)


def load_nearest_neighbors_model(chunk_size=1e6, fpformat="fp4", solubility='all', session=default_session,
                                 view=SequentialView(), model_size=100000, source="db", costum_query=None):
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

    if source == "file":
        load_nearest_neighbors_model_from_file(chunk_size=chunk_size, fpformat=fpformat, solubility=solubility,
                                               view=view, model_size=model_size)
    else:
        load_nearest_neighbors_model_from_db(fpformat=fpformat, solubility=solubility,
                                             model_size=model_size, session=session, costum_query=costum_query)


def load_nearest_neighbors_model_from_db(fpformat="fp4", solubility='all', model_size=1000, session=default_session,
                                         custom_query=None):
    """
    Loads a NN model.

    If a 'default_model.pickle' exists in data it will load the model. Otherwise it will build a model from the
    Database. This can take several hours depending on the size of the database.

    Parameters
    ----------
    fpformat : str
        The format of the fingerprint (see pybel.fps)
    solubility : str
        One of high, medium, low or all.
    model_size : int
        The size of each NearestNeighbor in the ensemble.
    session : Session
        SQLAlchemy session.
    custom_query : ClauseElement
        A query to filter elements from the database.

    """

    if custom_query is not None:
        indices = np.array(session.query(Metabolite.inchi_key).filter(custom_query).all(), dtype=INCHI_KEY_TYPE)
    else:
        indices = np.array(session.query(Metabolite.inchi_key).all(), dtype=INCHI_KEY_TYPE)

    n_models = math.ceil(len(indices) / model_size)
    chunk_size = math.ceil(len(indices) / n_models)
    chunks = [((i - 1) * chunk_size, i * chunk_size) for i in range(1, n_models + 1)]
    models = []

    for start, end in chunks:
        models.append(DBNearestNeighbors(indices[start:end], session, fpformat))

    return DistributedNearestNeighbors(models)


def load_nearest_neighbors_model_from_file(chunk_size=1e6, fpformat="fp4", solubility='all',
                                           view=SequentialView(), model_size=100000):
    """
    Loads a NN model from file.

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
        _indices, _features, _lengths = build_feature_table(Database.metabolites,
                                                            chunk_size=chunk_size,
                                                            fpformat=fpformat,
                                                            solubility=solubility,
                                                            view=view)
        pickle_large((_indices, _features, _lengths), model_file, progress=True)

    n_models = math.ceil(len(_indices) / model_size)
    nn_model = _build_nearest_neighbors_model(_indices, _features, _lengths, n_models)
    return nn_model


class DataBuilder(multiprocessing.Process):
    def __init__(self, inchi, similarity_cut, task_queue, results_queue, atoms_weight, bonds_weight, timeout,
                 *args, **kwargs):
        super(DataBuilder, self).__init__(*args, **kwargs)
        self._inchi = inchi
        self._similarity_cut = similarity_cut
        self._session = None
        self._tasks = task_queue
        self._results = results_queue
        self._atoms_weight = atoms_weight
        self._bonds_weight = bonds_weight
        self._timeout = timeout

    @property
    def session(self):
        if self._session is None:
            engine = create_engine('postgresql://', pool=connection_pool)
            session_maker = sessionmaker(engine)
            self._session = session_maker()
        return self._session

    @property
    def molecule(self):
        return rdkit.inchi_to_molecule(self._inchi)

    def run(self):
        while True:
            try:
                inchi_key, distance = self._tasks.get()
            except Empty:
                break
            else:
                result = self.apply_similarity_cut(inchi_key, distance)
                self._results.put(result)

        if self._session is not None:
            self.session.close()
            self.session.bind.dispose()

    def apply_similarity_cut(self, inchi_key, distance):
        met = Metabolite.get(inchi_key=inchi_key, session=self.session)

        try:
            molecule = met.molecule('rdkit', get3d=False)
            structural_similarity = rdkit.structural_similarity(self.molecule, molecule,
                                                                atoms_weight=self._atoms_weight,
                                                                bonds_weight=self._bonds_weight,
                                                                timeout=self._timeout)
            if structural_similarity >= self._similarity_cut:
                return [inchi_key, met.formula, met.num_atoms, met.num_bonds, 1 - distance, structural_similarity]
            else:
                return None
        except Exception as e:
            print(e)
            return None


def search_closest_compounds(molecule, nn_model=None, fp_cut=0.5, fpformat="maccs", atoms_diff=3,
                             bonds_diff=3, rings_diff=2, similarity_cut=0.6, session=default_session,
                             atoms_weight=0.5, bonds_weight=0.5, timeout=120):
    """
    Finds the closest compounds given a Molecule.


    Parameters
    ----------
    molecule : marsi.chemistry.molecule.Molecule
        A molecule representation.
    nn_model : marsi.nearest_neighbors.model.DistributedNearestNeighbors
        A nearest neighbors model.
    fp_cut : float
        A cutoff value for fingerprint similarity.
    fpformat : str
        A valid fingerprint format.
    atoms_diff : int
        The max number of atoms that can be different (in number, not type).
    bonds_diff : int
        The max number of bonds that can be different (in number, not type).
    rings_diff : int
        The max number of rings that can be different (in number, not type).
    similarity_cut : float
        A cutoff value for fingerprint similarity.
    session : Session
        SQLAlchemy session.
    atoms_weight : float
        The weight of having matching atoms in the structural similarity
    bonds_weight : float
        The weight of having matching bonds in the structural similarity

    Returns
    -------
    pandas.DataFrame
        A data frame with the closest InChI Keys as index and the properties calculated for each hit.
    """
    assert isinstance(molecule, Molecule)

    if nn_model is None:
        query = and_(Metabolite.num_atoms >= molecule.num_atoms - atoms_diff,
                     Metabolite.num_atoms <= molecule.num_atoms + atoms_diff,
                     Metabolite.num_bonds >= molecule.num_bonds - bonds_diff,
                     Metabolite.num_bonds <= molecule.num_bonds + bonds_diff,
                     Metabolite.num_rings >= molecule.num_rings - rings_diff,
                     Metabolite.num_rings <= molecule.num_rings + rings_diff)

        nn_model = load_nearest_neighbors_model_from_db(fpformat=fpformat, custom_query=query, session=session)

    assert isinstance(nn_model, DistributedNearestNeighbors)

    tasks_queue = multiprocessing.Queue()
    results_queue = multiprocessing.Queue()

    results = []

    neighbors = nn_model.radius_nearest_neighbors(molecule.fingerprint(fpformat), radius=fp_cut)

    dataframe = DataFrame(columns=["formula", "atoms", "bonds", "tanimoto_similarity", "structural_similarity"])

    progress = ProgressBar(maxval=len(neighbors), widgets=["Processing Neighbors: ", Bar(), ETA()])

    jobs = []
    for i in range(multiprocessing.cpu_count()):
        job = DataBuilder(molecule.inchi, similarity_cut, tasks_queue, results_queue,
                          atoms_weight, bonds_weight, timeout)
        jobs.append(job)
        job.start()

    for inchi_key, distance in six.iteritems(neighbors):
        tasks_queue.put((inchi_key, distance))

    progress.start()
    print(len(neighbors))
    while len(results) < len(neighbors):
        try:
            res = results_queue.get(timeout=500)
        except Empty:
            continue
        else:
            results.append(res)
            if res is not None:
                dataframe.loc[res[0]] = res[1:]
            progress.update(len(results))

    progress.finish()

    for job in jobs:
        job.terminate()

    return dataframe
