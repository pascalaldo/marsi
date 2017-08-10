# Copyright 2017 Chr. Hansen A/S and The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from cameo.core.target import Target
from cameo.flux_analysis.analysis import find_essential_metabolites
from gnomic import genotype_to_string, Genotype, Feature, Type, Ins

from marsi.cobra.flux_analysis.manipulation import knockout_metabolite, apply_anti_metabolite
from marsi.utils import search_metabolites


class AntiMetaboliteManipulationTarget(Target):
    """
    Metabolite target used to manipulate fluxes.
    There are three modes of action:
    1. Inhibition
    2. Competition
    3. Knockout (see Metabolite Knockout Target)

    If the metabolite is essential, the fluxes around the metabolite will be increased for the cells to survive. If the
    metabolite is not essential, the fluxes will be decreased as a consequence of the metabolite presence.

    """

    __gnomic_feature_type__ = "antimetabolite"

    def __init__(self, species_id, fraction=0.5, ignore_transport=True, allow_accumulation=True):
        super(AntiMetaboliteManipulationTarget, self).__init__(species_id)
        self.fraction = fraction
        self.ignore_transport = ignore_transport
        self.allow_accumulation = allow_accumulation

    def get_model_target(self, model):
        """
        Finds the metabolites in the model across compartments.

        Parameters
        ----------
        model: cameo.SolverBasedModel
            A constraint-based model.
        Returns
        -------
        list
            A list of cobra.core.metabolite.Metabolite.
        """
        return search_metabolites(model, self.id)

    def apply(self, model, reference=None):
        essential_metabolites = find_essential_metabolites(model)
        target_metabolites = self.get_model_target(model)

        apply_anti_metabolite(model, target_metabolites, essential_metabolites, reference,
                              allow_accumulation=self.allow_accumulation,
                              competition_fraction=self.fraction,
                              inhibition_fraction=self.fraction)

    def _repr_html_(self):  # pragma: no cover
        return "&#x2623;(%.3f)-%s" % (self.fraction, self.id)

    def __str__(self):
        return genotype_to_string(Genotype([self.to_gnomic()]))

    def __repr__(self):
        return "<AntiMetaboliteManipulation %s (%.3f)>" % (self.id, self.fraction)

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        feature = Feature(name=self.id, accession=accession, type=Type(self.__gnomic_feature_type__),
                          variant="value=%.5f" % self.fraction)
        return Ins(feature)


class MetaboliteKnockoutTarget(AntiMetaboliteManipulationTarget):
    def __init__(self, species_id, ignore_transport=True, allow_accumulation=True):
        super(MetaboliteKnockoutTarget, self).__init__(species_id, 1.0, ignore_transport, allow_accumulation)

    def apply(self, model, reference=None):
        for metabolite in self.get_model_target(model):
            knockout_metabolite(model, metabolite, self.ignore_transport, self.allow_accumulation)

    def _repr_html_(self):  # pragma: no cover
        return "&#x2623;%s" % self.id

    def __repr__(self):
        return "<MetaboliteKnockout %s>" % self.id

    def to_gnomic(self):
        accession = Target.to_gnomic(self)
        feature = Feature(name=self.id, accession=accession, type=Type(self.__gnomic_feature_type__))
        return Ins(feature)
