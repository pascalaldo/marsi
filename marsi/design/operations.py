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

from marsi.io.mongodb import Database
from marsi.processing.chemistry.openbabel import inchi_to_molecule, inchi_to_inchi_key

__all__ = ['AntiMetaboliteFilter', 'CombinedSubstratesAntiMetaboliteFilter']


class AntiMetaboliteFilter(object):
    """
    Filter's antimetabolites.

    The filter consists of multiple steps of filtering:
        1. Fingerprint similarity of radios _min_tanimoto_
        2. Solubility > min_solubility



    """
    def __init__(self, nn_model, metabolites=Database.metabolites.collection):
        self.nn_model = nn_model
        self.metabolites = metabolites

    def __call__(self, metabolite, min_tanimoto=0.5, min_solubility=0.00006, fpformat='ecfp10', mode="cl"):
        """
        Arguments
        ---------

        metabolite: cobra.core.Metabolite
            A metabolite to search.
        min_tanimoto: float
            Tanimoto similarity cutoff.
        min_solubility: float
            Solubility cutoff.
        fpformat: str
            Fingerprint format (one of pybel.fps).
        mode: str
            'cl' to run the openCL implementation of the fingerprint search.

        Returns
        -------
        tuple
            InChI, InChI Key, list of matches.
        """
        from marsi.processing.models import find_inchi_for_bigg_metabolite
        try:
            try:
                inchi = metabolite.annotation['inchi']
            except KeyError:
                metabolite.annotation['inchi'] = inchi = find_inchi_for_bigg_metabolite(metabolite.model, metabolite.id)

            inchi_key = inchi_to_inchi_key(inchi)

            filter_1 = self._filter_by_fingerprint(inchi, inchi_key, fpformat, min_tanimoto, mode)
            filter_2 = self._filter_native_metabolites(filter_1.keys(), metabolite.model)
            filter_3 = self._filter_by_solubility(filter_2, min_solubility)
            return inchi, inchi_key, list(filter_3.keys())
        except ValueError:
            return None, None, []

    def _filter_by_fingerprint(self, inchi, inchi_key, fpformat, min_tanimoto, mode):
        fingerprint = self._calculate_fingerprint(inchi, inchi_key, fpformat)
        return self.nn_model.radius_nearest_neighbors(fingerprint, radius=1-min_tanimoto, mode=mode)

    def _calculate_fingerprint(self, inchi, inchi_key, fpformat):
        try:
            metabolite = self.metabolites.get(inchi_key)
            fingerprint = metabolite.fingerprint(fpformat)
        except KeyError:
            metabolite = inchi_to_molecule(inchi)
            fingerprint = metabolite.calcfp(fpformat).fp

        return fingerprint

    @staticmethod
    def _filter_native_metabolites(inchi_keys, model):
        model_inchi = [m.annotation.get('inchi', None) for m in model.metabolites]
        model_inchi_keys = [inchi_to_inchi_key(inchi) for inchi in model_inchi if inchi]
        return [inchi_key for inchi_key in inchi_keys if inchi_key in model_inchi_keys]

    def _filter_by_solubility(self, inchi_keys, min_solubility):
        solubility = {key: self.metabolites.get(key).solubility for key in inchi_keys}
        return {key: sol for key, sol in solubility.items() if sol >= min_solubility}


class CombinedSubstratesAntiMetaboliteFilter(AntiMetaboliteFilter):
    def __call__(self, metabolites, min_tanimoto=0.5, min_solubility=5, fpformat='ecfp4', mode="cl"):
        from marsi.processing.models import find_inchi_for_bigg_metabolite
        try:
            inchis = [find_inchi_for_bigg_metabolite(m.model, m.id) for m in metabolites]
            inchi_keys = [inchi_to_inchi_key(inchi) for inchi in inchis]
            return inchis, inchi_keys
        except ValueError:
            return None, None, []


class ReplaceKnockoutFinder(object):
    def __init__(self, model, simulation_method, simulation_kwargs, currency_metabolites):
        self.model = model
        self.simulation_method = simulation_method
        self.simulation_kwargs = simulation_kwargs
        self.currency_metabolites = currency_metabolites

    def __call__(self, targets, fitness, objective_function, cache, max_lost=0.2):
        from marsi.processing.models import test_reaction_knockout_replacements
        simulation_kwargs = dict(cache=cache)
        simulation_kwargs.update(self.simulation_kwargs)
        return test_reaction_knockout_replacements(self.model, targets, objective_function, fitness,
                                                   self.simulation_method, simulation_kwargs=simulation_kwargs,
                                                   ignore_metabolites=self.currency_metabolites, max_lost=max_lost)
