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

ctypedef np.int32_t INT32_t
ctypedef np.float32_t FLOAT32_t

cdef float _tanimoto_coefficient(np.ndarray[INT32_t, ndim=1] fingerprint1, np.ndarray[INT32_t, ndim=1] fingerprint2)
cdef float _tanimoto_distance(np.ndarray[INT32_t, ndim=1] fingerprint1, np.ndarray[INT32_t, ndim=1] fingerprint2)
