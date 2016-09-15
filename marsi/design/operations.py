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
from numpy import array

from marsi.distances import tanimoto_distance
from marsi.processing.models import find_bigg_metabolite
from marsi.processing.chemistry import inchi_to_molecule

__all__ = ['NearestNeighbors']


class NearestNeighbors(object):
    def __init__(self, database):
        self.database = database

    def __call__(self, metabolite, min_tanimoto=0.8, min_rmsd=0.8, fpformat='maccs'):
        try:
            inchi, inchi_key = find_bigg_metabolite(metabolite.model, metabolite.id)
            metabolite_molecule = inchi_to_molecule(inchi)
            metabolite_fingerprint = metabolite_molecule.calcfp(fpformat)
            fingerprints = [m.fingerprint(fpformat=fpformat) for m in self.database.metabolites]
            tamimoto_distances = array([tanimoto_distance(metabolite_fingerprint, fp) for fp in fingerprints])
            tier1 = tamimoto_distances[tamimoto_distances < min_tanimoto]

            # rmsd_distances = tier1.apply(rsmd, args=(inchi, ))
            # tier2 = rmsd_distances[rmsd_distances < rmsd_distance]
            return inchi, inchi_key, tier1
        except ValueError:
            return None, None, []

