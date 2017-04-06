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

from numpy import ndarray

from libc.math cimport sqrt

cdef extern int __builtin_popcount(unsigned int) nogil

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.cdivision(True)
cdef float _tanimoto_coefficient(np.ndarray[INT32_t, ndim=1] fingerprint1, np.ndarray[INT32_t, ndim=1] fingerprint2):
    """
    Calculate the Tanimoto coefficient for 2 fingerprints.

    Parameters
    ----------
    fingerprint1 : ndarray
        First fingerprint.
    fingerprint2 : ndarray
        Second fingerprint.

    Returns
    -------
    float
        The Tanimoto coefficient.
    """
    cdef unsigned int len1 = fingerprint1.shape[0]
    cdef unsigned int len2 = fingerprint2.shape[0]

    if len1 != len2:
        return -1.0

    cdef int and_bits = 0
    cdef int or_bits = 0

    cdef int fp_and
    cdef int fp_or

    for i in range(len1):
        fp_and = fingerprint1[i] & fingerprint2[i]
        fp_or = fingerprint1[i] | fingerprint2[i]

        and_bits += __builtin_popcount(fp_and)
        or_bits += __builtin_popcount(fp_or)

    return float(and_bits)/float(or_bits)

@cython.nonecheck(False)
@cython.cdivision(True)
cdef float _tanimoto_distance(np.ndarray[INT32_t, ndim=1] fingerprint1, np.ndarray[INT32_t, ndim=1] fingerprint2):
    """
    Calculate the Tanimoto distance for 2 fingerprints (1 - tanimoto coefficient).

     Arguments
    ---------
    fingerprint1: ndarray
        First fingerprint.
    fingerprint2: ndarray
        Second fingerprint.

    Returns
    -------
    float
        The Tanimoto distance.
    """

    return 1 - _tanimoto_coefficient(fingerprint1, fingerprint2)


def tanimoto_coefficient(np.ndarray[INT32_t, ndim=1] fingerprint1, np.ndarray[INT32_t, ndim=1] fingerprint2):
    """
    Calculate the Tanimoto coefficient for 2 fingerprints.

    Parameters
    ----------
    fingerprint1 : ndarray
        First fingerprint.
    fingerprint2 : ndarray
        Second fingerprint.

    Returns
    -------
    float
        The Tanimoto coefficient.
    """
    return _tanimoto_coefficient(fingerprint1, fingerprint2)


def tanimoto_distance(np.ndarray[INT32_t, ndim=1] fingerprint1, np.ndarray[INT32_t, ndim=1] fingerprint2):
    """
    Calculate the Tanimoto distance for 2 fingerprints (1 - tanimoto coefficient).

    Parameters
    ----------
    fingerprint1 : ndarray
        First fingerprint.
    fingerprint2 : ndarray
        Second fingerprint.

    Returns
    -------
    float
        The Tanimoto distance.
    """

    return _tanimoto_distance(fingerprint1, fingerprint2)


def rmsd(np.ndarray[FLOAT32_t, ndim=3] v, np.ndarray[FLOAT32_t, ndim=3] w):
    """
    Root-mean-squared deviation of XYZ.

    $$RMSD = \sqrt{\frac{1}{n} \sum_{i=1}^{n}{(v - t)^2}}}$$

    Parameters
    ----------
    v : list
        List of x, y, z
    w : list
        List of x, y, z
    """

    if len(v) != len(w):
        raise ValueError("Length of vectors should be the same")

    cdef unsigned int n = len(v)
    cdef unsigned int sum = 0

    for i in range(n):
        sum += (v[i] - w[i])**2 + (v[i] - w[i])**2 + (v[i] - w[i])**2

    return sqrt(1/n * sum)

@cython.nonecheck(False)
@cython.cdivision(True)
cdef limits(np.ndarray[FLOAT32_t, ndim=2] coords, np.ndarray[FLOAT32_t, ndim=1] vdw_radii, int index):
    cdef unsigned int i = 0
    cdef unsigned int n = len(coords)
    cdef np.ndarray[FLOAT32_t, ndim=1] _coord = np.zeros(n, dtype=np.float32)

    for i in range(len(coords)):
        _coord[i] = coords[i][index] + vdw_radii[i]
    return min(_coord), max(_coord)

@cython.nonecheck(False)
@cython.cdivision(True)
cdef unsigned int point_in_molecule(np.ndarray[FLOAT32_t, ndim=2] coords, np.ndarray[FLOAT32_t, ndim=1] vdw_radii, np.ndarray[FLOAT32_t, ndim=1] point):

    cdef float r, x, y, z
    cdef unsigned int n = len(coords)
    cdef np.ndarray[FLOAT32_t, ndim=1] _coord
    for i in range(n):
        _coord = coords[i]
        r = vdw_radii[i]

        x = point[0] - _coord[0]
        y = point[1] - _coord[1]
        z = point[2] - _coord[2]

        if x*x + y*y + z*z <= r**2:
            return 1
    return 0

@cython.nonecheck(False)
@cython.cdivision(True)
cdef _monte_carlo_volume(np.ndarray[FLOAT32_t, ndim=2] coords, np.ndarray[FLOAT32_t, ndim=1] vdw_radii, float tolerance,
                         unsigned int max_iterations, unsigned int step_size, int seed, unsigned int verbose):

    np.random.seed(seed)

    # determine the bounding box
    cdef float[2] mol_x = limits(coords, vdw_radii, 0)
    cdef float[2] mol_y = limits(coords, vdw_radii, 2)
    cdef float[2] mol_z = limits(coords, vdw_radii, 2)

    cdef float mol_x_len = mol_x[1] - mol_x[0]
    cdef float mol_y_len = mol_y[1] - mol_y[0]
    cdef float mol_z_len = mol_y[1] - mol_y[0]

    cdef float mol_x_mid = (mol_x[0] + mol_x[1]) / 2.0
    cdef float mol_y_mid = (mol_y[0] + mol_y[1]) / 2.0
    cdef float mol_z_mid = (mol_z[0] + mol_z[1]) / 2.0

    cdef float box_x_min = mol_x_mid - 2.0 * mol_x_len
    cdef float box_x_max = mol_x_mid + 2.0 * mol_x_len
    cdef float box_y_min = mol_y_mid - 2.0 * mol_y_len
    cdef float box_y_max = mol_y_mid + 2.0 * mol_y_len
    cdef float box_z_min = mol_z_mid - 2.0 * mol_z_len
    cdef float box_z_max = mol_z_mid + 2.0 * mol_z_len

    cdef float box_x_len = box_x_max - box_x_min
    cdef float box_y_len = box_y_max - box_y_min
    cdef float box_z_len = box_z_max - box_z_min
    cdef float box_volume = box_x_len * box_y_len * box_z_len

    if verbose:
        print("Box volume %.5f" % box_volume)

    # calculate the ratio of points inside and outside the VdW shell of the molecule

    # one hundred thousand points to start
    cdef unsigned int total_points = 100000
    cdef unsigned int points_in_molecule = 0

    cdef np.ndarray[FLOAT32_t, ndim=2] points = np.random.uniform(box_x_min, box_x_max, (total_points, 3)).astype(np.float32)
    cdef np.ndarray[FLOAT32_t, ndim=1] point
    cdef int i =0
    for i in range(total_points):
        point = points[i]
        points_in_molecule += point_in_molecule(coords, vdw_radii, point)

    # Adding more points until reach some convergence
    cdef float volume = box_volume * float(points_in_molecule) / float(total_points)
    cdef float new_volume = 0.0
    cdef j = 0
    while abs(volume - new_volume) > tolerance and j <= max_iterations:
        points = np.random.uniform(box_x_min, box_x_max, (step_size, 3)).astype(np.float32)
        for i in range(step_size):
            point = points[i]
            points_in_molecule += point_in_molecule(coords, vdw_radii, point)
        total_points += step_size
        volume = new_volume
        new_volume =  box_volume * float(points_in_molecule) / float(total_points)
        if verbose:
            print("Iteration %i, volume: %.5f, new volume: %.5f" % (j, volume, new_volume))
        j += 1

    return volume

def monte_carlo_volume(np.ndarray[FLOAT32_t, ndim=2] coords, np.ndarray[FLOAT32_t, ndim=1] vdw_radii, float tolerance,
                         int max_iterations, int step_size, int seed, int verbose):
    """
    Adapted from:

    Simple Monte Carlo estimation of VdW molecular volume (in A^3)
    by Geoffrey Hutchison <geoffh@pitt.edu>

    https://github.com/ghutchis/hutchison-cluster
    """

    return _monte_carlo_volume(coords, vdw_radii, tolerance, max_iterations, step_size, seed, verbose)
