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


__all__ = ['NearestNeighbors']


class NearestNeighbors(object):
    def __init__(self, database):
        self.database = database

    def __call__(self, metabolite, tanimoto_distance=0.8, rmsd_distance=0.8):
        inchi, inchi_key = find_model_metabolite(metabolite.id)

        tamimoto_distances = self.database.metabolites.inchi.apply(tanimoto, args=(inchi, ))
        tier1 = tamimoto_distances[tamimoto_distances.distance < tanimoto_distance]

        rmsd_distances = tier1.apply(rsmd, args=(inchi, ))
        tier2 = rmsd_distances[rmsd_distances < rmsd_distance]
        return inchi, inchi_key, tier2
