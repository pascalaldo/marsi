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

from cameo.parallel import SequentialView

try:  # pragma: no cover
    import pyopencl as cl
    MF = cl.mem_flags
    cl_available = True
except ImportError:
    class ClFail:
        def __getattr__(self, item):
            raise RuntimeError("OpenCL is not available. try 'pip install pyopencl'")

    cl = ClFail()
    MF = None
    cl_available = False

import numpy as np
from pandas import DataFrame

from marsi.utils import timing
from marsi.nearest_neighbors import model_ext

# Tanimoto coefficient calculation is implemented based on to OpenBabel's implementation
# https://github.com/openbabel/openbabel/blob/master/include/openbabel/fingerprint.h#L86
nn_source = """
__kernel void nn(__global float *distances, __global int *fingerprint, __constant int *fpLength,
                 __global int *database, __global int *positions, __constant int *lengths) {
    int gid = (int)get_global_id(0);

    float tanimotoCoefficient = -1;

    int dbEntryLength = lengths[gid];
    int fpLen = fpLength[0];

    if (fpLen == dbEntryLength) {
        int andBits = 0;
        int orBits = 0;
        int fpPos = positions[gid];
        int i;
        for (i=0; i < fpLen; ++i) {
            int fpAnd = database[fpPos+i] & fingerprint[i];
            int fpOr = database[fpPos+i] | fingerprint[i];
            andBits += popcount(fpAnd);
            orBits += popcount(fpOr);
        }

        tanimotoCoefficient = (float)andBits/(float)orBits;
    }
    distances[gid] = 1 - tanimotoCoefficient;
}
"""


class KNN(object):
    """
    K-Nearest Neighbors runner object.

    It is assigned to a model and runs the `knn` function.

    Attributes
    ----------
    fp : numpy.array
        A numpy.array with the fingerprint values.
    k : int
        The maximum number of neighbors to retrieve.
    mode : str
        'native' to run python implementation or 'cl' to run OpenCL implementation if available.
    """
    def __init__(self, fingerprint, k, mode):
        self.fp = np.array(list(fingerprint), dtype=np.int32)
        self.k = k
        self.mode = mode

    def __call__(self, nn):
        return nn.knn(self.fp, k=self.k)

    def __getstate__(self):
        return dict(fp=self.fp.tolist(), k=self.k, mode=self.mode)

    def __setstate__(self, d):
        d['fp'] = np.array(d['fp'], dtype=np.int32)
        self.__dict__.update(d)


class RNN(object):
    """
    R-Nearest Neighbors runner object.

    It is assigned to a model and runs the `rnn` function.

    Attributes
    ----------
    fp : numpy.array
        A numpy.array with the fingerprint values.
    radius : float
        A distance radius ]0, 1].
    mode : str
        'native' to run python implementation or 'cl' to run OpenCL implementation if available.
    """
    def __init__(self, fingerprint, radius, mode):
        self.fp = np.array(list(fingerprint), dtype=np.int32)
        assert 0 < radius <= 1
        self.radius = radius
        self.mode = mode

    def __call__(self, nn):
        return nn.rnn(self.fp, radius=self.radius, mode=self.mode)

    def __getstate__(self):
        return dict(fp=self.fp.tolist(), radius=self.radius, mode=self.mode)

    def __setstate__(self, d):
        d['fp'] = np.array(d['fp'], dtype=np.int32)
        self.__dict__.update(d)


class Distance(object):
    """
    Distance runner object.

    It is assigned to a model and runs the `distance` function.

    Attributes
    ----------
    fp : numpy.array
        A numpy.array with the fingerprint values.
    mode : str
        'native' to run python implementation or 'cl' to run OpenCL implementation if available.
    """
    def __init__(self, fingerprint, mode):
        self.fp = np.array(list(fingerprint), dtype=np.int32)
        self.mode = mode

    def __call__(self, nn):
        return nn.distances(self.fp, mode=self.mode)

    def __getstate__(self):
        return dict(fp=self.fp, mode=self.mode)

    def __setstate__(self, d):
        d['fp'] = np.array(d['fp'], dtype=np.int32)
        self.__dict__.update(d)


class DistributedNearestNeighbors(object):
    """
    Nearest Neighbors distributed implementation.

    Attributes
    ----------
    index : numpy.array
        The index of all entries across multiple models.
    """
    def __init__(self, nns):
        self._nns = nns

    def k_nearest_neighbors(self, fingerprint, k=5, mode="native", view=SequentialView()):
        """
        Retrieves the K nearest neighbors to a fingerprint.

        Parameters
        ----------
        fingerprint : list, np.array, tuple
            A fingerprint to use as query.
        k : int
            The number of neighbors to retrieve.
        mode : str
            'native' to run python implementation or 'cl' to run OpenCL implementation if available.
        view : cameo.parallel.ParallelView, cameo.parallel.SequentialView
            A parallel mode runner.
        Returns
        -------
        dict
            A dictionary with the InChI Key as key and the distance as value.
        """
        func = KNN(fingerprint, k, mode)
        results = view.map(func, self._nns)
        neighbors = {}
        [neighbors.update(res) for res in results]
        return dict(sorted(neighbors.items(), key=lambda x: x[1])[:k])

    def radius_nearest_neighbors(self, fingerprint, radius=0.25, mode="native", view=SequentialView()):
        """
        Retrieves the nearest neighbors to a fingerprint within a distance radius.

        Parameters
        ----------
        fingerprint : list, np.array, tuple
            A fingerprint to use as query.
        radius : float
            A distance radius ]0, 1].
        mode : str
            'native' to run python implementation or 'cl' to run OpenCL implementation if available.
        view : cameo.parallel.ParallelView, cameo.parallel.SequentialView
            A parallel mode runner.
        Returns
        -------
        dict
            A dictionary with the InChI Key as key and the distance as value.
        """
        func = RNN(fingerprint, radius, mode)
        results = view.map(func, self._nns)
        neighbors = {}
        [neighbors.update(res) for res in results]
        return neighbors

    def distances(self, fingerprint, mode="native", view=SequentialView()):
        """
        Retrieves the distance a fingerprint and all elements in the model.

        Parameters
        ----------
        fingerprint : list, np.array, tuple
            A fingerprint to use as query.
        mode : str
            'native' to run python implementation or 'cl' to run OpenCL implementation if available.
        view : cameo.parallel.ParallelView, cameo.parallel.SequentialView
            A parallel mode runner.
        Returns
        -------
        dict
            A dictionary with the InChI Key as key and the distance as value.
        """
        func = Distance(fingerprint, mode)
        results = view.map(func, self._nns)
        distances = {}
        [distances.update(res) for res in results]
        return distances

    def __repr__(self):
        return " | ".join(repr(nn) for nn in self._nns)

    def __getitem__(self, index):
        return self._nns[index]

    @property
    def index(self):
        return np.concatenate([nn.index for nn in self._nns])

    def distance_matrix(self, mode="native"):
        """
        Generates a distance matrix between all elements in the models.

        Parameters
        ----------
        mode : str
            'native' to run python implementation or 'cl' to run OpenCL implementation if available.

        Returns
        -------
        numpy.array
            The distance matrix.
        """
        index = self.index
        size = len(index)
        matrix = np.zeros((size, size), dtype=np.float)
        for i in range(size):
            matrix[i] = self.distances(self.feature(i), mode=mode)

        return matrix

    def feature(self, index):
        """
        Retrieves the fingerprint at a given index. The index is global for the ensemble of models.

        Returns
        -------
        numpy.array
            The fingerprint.

        Raises
        ------
        IndexError
        """
        total = len(self)
        if index >= total:
            raise IndexError(index)

        group = 0
        off = 0
        for nn in self._nns:
            if index + 1 > len(nn):
                group += 1
                off += len(nn)

        return self._nns[group][index - off]

    def __len__(self):
        return sum(len(nn) for nn in self._nns)


class NearestNeighbors(model_ext.CNearestNeighbors):
    def __init__(self, index, features, features_lengths, use_cl=False, opencl_context=None):
        features = np.concatenate(features).astype(np.int32)
        features_lengths = np.array(features_lengths, dtype=np.int32)
        super(NearestNeighbors, self).__init__(features, features_lengths)
        self._index = index
        self._use_cl = use_cl
        if cl_available and self._use_cl:
            if opencl_context is None:
                self._ctx = cl.create_some_context()
            else:
                self._ctx = opencl_context
            self._program = cl.Program(self._ctx, nn_source).build()

    @property
    def cl_context(self):
        if not self._use_cl:
            raise RuntimeError("Configuration is set to don't run CL")
        return self._ctx

    @cl_context.setter
    def cl_context(self, ctx):
        if not self._use_cl:
            raise RuntimeError("Configuration is set to don't run CL")
        self._ctx = ctx
        self._program = cl.Program(self._ctx, nn_source).build()

    def queue(self):
        return cl.CommandQueue(self.cl_context)

    @property
    def program(self):
        if not self._use_cl:
            raise RuntimeError("Configuration is set to don't run CL")
        return self._program

    def __getstate__(self):
        state = super(NearestNeighbors, self).__getstate__()
        state.update({"_index": self._index, "_use_cl": self._use_cl})
        return state

    def __getitem__(self, index):
        start = self.start_positions[index]
        length = self.start_positions[index]
        end = start + length
        return self._features[start:end]

    def __setstate__(self, state):
        super(NearestNeighbors, self).__setstate__(state)
        self._index = state['_index']
        self._use_cl = state['_use_cl']
        if cl_available and self._use_cl:
            self.cl_context = cl.create_some_context()

    def knn(self, fingerprint, k, mode="native"):
        """
        K-Nearest Neighbors

        Parameters
        ----------

        fingerprint : ndarray
            The fingerprint to search for.
        k : int
            The number of neighbors to return.
        mode : str
            'native' to run python implementation or 'cl' to run OpenCL implementation if available.

        Returns
        -------
        dict
            (Index --> Distance)

        """
        distances = self.distances(fingerprint, mode)
        indices = np.argsort(distances)
        return {i.decode('utf-8'): d for i, d in zip(self._index[indices, 0][:k], distances[indices][:k])}

    def rnn(self, fingerprint, radius, mode="native"):
        """
        Radius-Nearest Neighbors

        Parameters
        ----------
        fingerprint : ndarray
            The fingerprint to search for.
        radius : float
            The maximum distance of neighbors to return.
        mode : str
            'native' to run python implementation or 'cl' to run OpenCL implementation if available.

        Returns
        -------
        dict
            (Index --> Distance)

        """
        distances = self.distances(fingerprint, mode)
        indices = np.argsort(distances)
        indices = indices[distances[indices] <= radius]
        return {i.decode('utf-8'): d for i, d in zip(self._index[indices, 0], distances[indices])}

    def distances(self, fingerprint, mode="native"):
        if mode == "native":
            return self.distances_py(fingerprint)
        elif mode == "cl":
            return self.distances_cl(fingerprint)
        else:
            raise ValueError("'mode' can only be 'native' or 'cl'")

    def max_memory_allocation_size(self):
        return min([d.max_mem_alloc_size for d in self.cl_context.devices])

    @timing(debug=True)
    def distances_cl(self, fingerprint):
        database_buffer = self.input_buffer(self.features)

        start_positions_buffer = self.input_buffer(self.start_positions)

        features_lengths_buffers = self.input_buffer(self.features_lengths)

        fingerprint_buffer = self.input_buffer(fingerprint)

        fingerprint_length = np.array([len(fingerprint)])
        fingerprint_length_buffer = self.input_buffer(fingerprint_length)

        distances = np.zeros(len(self._index), dtype=np.float32)
        distances_buf = self.output_buffer(distances)

        # Run the kernel
        queue = self.queue()
        exec_evt = self.run_kernel(queue, distances.shape, distances_buf, fingerprint_buffer, fingerprint_length_buffer,
                                   database_buffer, start_positions_buffer, features_lengths_buffers)

        # Read the result
        cl.enqueue_copy(queue, distances, distances_buf, is_blocking=True, wait_for=[exec_evt]).wait()
        queue.finish()
        return distances

    @timing(debug=True)
    def input_buffer(self, array):
        return cl.Buffer(self._ctx, MF.READ_ONLY | MF.COPY_HOST_PTR, hostbuf=array)

    @timing(debug=True)
    def output_buffer(self, array):
        return cl.Buffer(self._ctx, MF.WRITE_ONLY, size=array.nbytes)

    @timing(debug=True)
    def run_kernel(self, queue, shape, *args):
        exec_evt = self.program.nn(queue, shape, None, *args)
        return exec_evt

    @property
    def index(self):
        return self._index

    @property
    def data_frame(self):
        """
        Create a DataFrame with this model.

        Returns
        -------
        pandas.DataFrame
            A data frame.
        """
        df = DataFrame([self[i] for i in range(len(self))],
                       index=["fingerprint"], columns=self._index)
        return df.T

    def __repr__(self):
        return "NearestNeighbors %i rows" % (len(self))

    def __len__(self):
        return len(self._index)
