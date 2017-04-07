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
import numpy as np
from cachetools import LRUCache
from scipy.spatial import ConvexHull

from marsi.chemistry.common_ext import tanimoto_coefficient, tanimoto_distance, rmsd, monte_carlo_volume
from scipy.spatial.qhull import QhullError


__all__ = ["rmsd", "tanimoto_coefficient", "tanimoto_distance", "monte_carlo_volume",
           "INCHI_KEY_REGEX", 'SOLUBILITY']


inchi_key_lru_cache = LRUCache(maxsize=512)


SOLUBILITY = {
    "high": lambda sol: sol > 0.00006,
    "medium": lambda sol: 0.00001 <= sol <= 0.00006,
    "low": lambda sol: sol < 0.00001,
    "all": lambda sol: True
}


INCHI_KEY_REGEX = re.compile("[0-9A-Z]{14}\-[0-9A-Z]{8,10}\-[0-9A-Z]")


def convex_hull_volume(xyz):
    try:
        return ConvexHull(xyz).volume
    except (QhullError, ValueError):
        return np.nan