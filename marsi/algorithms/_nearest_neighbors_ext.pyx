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

import cython
cimport cython

import numpy as np
cimport numpy as np

from marsi.chemistry.common_ext cimport _tanimoto_distance

ctypedef np.int32_t INT32_t
ctypedef np.float32_t FLOAT32_t


cdef extern int __builtin_popcount(unsigned int) nogil

@cython.nonecheck(False)
@cython.cdivision(True)
cdef class CNearestNeighbors:
    cdef np.ndarray _features
    cdef np.ndarray _features_lengths
    cdef np.ndarray _start_positions

    def __init__(self, features, features_lengths):
        self._features = features
        self._features_lengths = features_lengths
        self._start_positions = None

    def __getstate__(self):
        return {
            '_features': self._features,
            '_features_lengths': self._features_lengths,
            '_start_positions': self._start_positions
        }

    def __setstate__(self, state):
        self._features = state['_features']
        self._features_lengths = state['_features_lengths']
        self._start_positions = state['_start_positions']

    @property
    def features(self):
        return self._features

    @property
    def features_lengths(self):
        return self._features_lengths

    @property
    def start_positions(self):
        if self._start_positions is None:
            self._start_positions = self._calc_start_positions()
        return self._start_positions

    def distances_py(self, fingerprint):
        return self._distances(fingerprint)

    cdef np.ndarray[INT32_t, ndim=1] _calc_start_positions(self):
        cdef unsigned int i = 0
        cdef unsigned int n = len(self._features_lengths)

        cdef np.ndarray[INT32_t, ndim=1] res = np.empty(n, dtype=np.int32)

        res[0] = i
        for pos in range(n-1):
            i += self._features_lengths[pos]
            res[pos+1] = i

        return res

    cdef np.ndarray[FLOAT32_t, ndim=1] _distances(self, np.ndarray[INT32_t, ndim=1] fingerprint):
        cdef unsigned int n = len(self._features_lengths)
        cdef np.ndarray[FLOAT32_t, ndim=1] _distances = np.zeros(n, dtype=np.float32)

        for i in range(n):
            start = self.start_positions[i]
            size = self._features_lengths[i]
            _distances[i] = _tanimoto_distance(fingerprint, self._features[start:start+size])

        return _distances

