# Copyright 2016 Chr. Hansen A/S
# Copyright 2016 The Novo Nordisk Foundation Center for Biosustainability, DTU.

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

# http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from cameo import fba
from cameo.strain_design import OptKnock
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield

from marsi.design.evolutionary import OptGeneMet
from marsi.design.operations import AntiMetaboliteFilter, CombinedSubstratesAntiMetaboliteFilter

CURRENCY_METABOLITES = ["atp", "adp", "nad", "nadh", "nadp", "nadph", "amp",
                        "h2o", "h", "coa", "acp", "pi", 'pppi', 'ppi']


class GenericMARSIDesignMethod(object):
    """
    Generic wrapper for Metabolite Analog design method.

    This one just runs a optimization method

    Example
    -------
    >>> from marsi import design
    >>> from cameo import models
    >>> from cameo.strain_design import OptGene
    >>> designer = design.GenericMARSIDesignMethod(model=models.bigg.iJO1366)
    >>> designer.optimize_with_reaction("succ_e", max_interventions=5, substrate="EX_glc__D_e",
    >>> biomass="BIOMASS_Ec_iJO1366_core_53p95M", max_results=25, design_method=OptGene, manipulation_type="reactions")


    """
    def __init__(self, model=None, nearest_neighbors_model=None, min_tanimoto=0.75, currency_metabolites=None):
        self.model = model
        self.nearest_neighbors_model = nearest_neighbors_model
        self.min_tanimoto = min_tanimoto
        self.currency_metabolites = currency_metabolites or CURRENCY_METABOLITES

    def map_metabolite_neighbors(self, metabolite_id, fpformat='fp4', mode='cl'):
        anti_metabolite_filter = AntiMetaboliteFilter(self.nearest_neighbors_model)
        metabolite = self.model.metabolites.get_by_id(metabolite_id)
        return anti_metabolite_filter(metabolite, min_tanimoto=self.min_tanimoto, fpformat=fpformat, mode=mode)

    def transition_state_neighbors(self, metabolite_ids, fpformat='fp4', mode='cl'):
        anti_metabolite_filter = CombinedSubstratesAntiMetaboliteFilter(self.nearest_neighbors_model)
        metabolites = [self.model.metabolites.get_by_id(mid) for mid in metabolite_ids]
        return anti_metabolite_filter(metabolites, min_tanimoto=self.min_tanimoto, fpformat=fpformat, mode=mode)

    def optimize_with_reaction(self, target, max_interventions=1, substrate=None,
                               biomass=None, design_method=OptKnock, max_results=100,
                               non_essential_metabolites=False, **design_kwargs):

        target_flux = self.model._reaction_for(target)

        exclude_reactions = []
        if non_essential_metabolites:
            exclude_reactions = self.essential_metabolites_reactions()

        if 'essential_reactions' in design_kwargs:
            design_kwargs['essential_reactions'] += exclude_reactions
        else:
            design_kwargs['essential_reactions'] = exclude_reactions

        bpcy = biomass_product_coupled_yield(biomass, target_flux, substrate)

        designer = design_method(model=self.model, **design_kwargs)
        knockouts = designer.run(max_knockouts=max_interventions, biomass=biomass, substrate=substrate,
                                 target=target_flux, max_results=max_results)

        from marsi.design.post_processing import replace_knockouts
        return replace_knockouts(self.model, knockouts.data_frame, bpcy, fba, {}, self.currency_metabolites)

    def optimize_with_metabolites(self, target, maximize=True, max_interventions=1, substrate=None, biomass=None,
                                  max_results=100, non_essential_metabolites=False, **design_kwargs):

        target_flux = self.model._reaction_for(target)
        designer = OptGeneMet(model=self.model, essential_metabolites=CURRENCY_METABOLITES, **design_kwargs)
        knockouts = designer.run(max_knockouts=max_interventions, biomass=biomass, substrate=substrate,
                                 target=target_flux, max_results=max_results)

        return knockouts.data_frame

    def essential_metabolites_reactions(self):
        essential_metabolites = self.model.essential_metabolites()
        reactions = set()
        for metabolite in essential_metabolites:
            reactions.update(metabolite.reactions)

        return reactions


class RandomMutagenesisDesign(GenericMARSIDesignMethod):
    """
    Apply only knockout like designs where total loss of functions are expected.



    """
    def __call__(self, target, maximize=True, **kwargs):
        pass


class ALEDesign(GenericMARSIDesignMethod):
    """
    Apply both knockout and flux modulation. The strains will be selected via growth rate proxy.


    """
    def __call__(self, target, maximize=True, **kwargs):
        pass
