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
from cameo.core import SolverBasedModel
from cameo.strain_design import DifferentialFVA
from cameo.strain_design import OptKnock, OptGene
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
from cameo.util import TimeMachine
from pandas import DataFrame

from marsi.cobra.strain_design.evolutionary import OptMet
from marsi.cobra.strain_design.post_processing import replace_design
from marsi.cobra.utils import CURRENCY_METABOLITES


class GenericMARSIDesignMethod(object):
    """
    Generic wrapper for Metabolite Analog design method.

    This one just runs a optimization method

    Example
    -------
    >>> from marsi.cobra import strain_design
    >>> from cameo import models
    >>> from cameo.strain_design import OptGene
    >>> designer = strain_design.GenericMARSIDesignMethod(model=models.bigg.iJO1366)
    >>> designer.optimize_with_reaction("succ_e", max_interventions=5, substrate="EX_glc__D_e",
    >>> biomass="BIOMASS_Ec_iJO1366_core_53p95M", max_results=25, design_method=OptGene, manipulation_type="reactions")


    """
    def __init__(self, model=None, nearest_neighbors_model=None, min_tanimoto=0.75, currency_metabolites=None,
                 essential_metabolites=None):
        assert isinstance(model, SolverBasedModel)
        self.model = model
        self.nearest_neighbors_model = nearest_neighbors_model
        self.min_tanimoto = min_tanimoto
        self.currency_metabolites = currency_metabolites or CURRENCY_METABOLITES
        if essential_metabolites is None:
            essential_metabolites = []

        self.essential_metabolites = model.essential_metabolites() + essential_metabolites
        for i, m in enumerate(self.essential_metabolites):
            if isinstance(m, str):
                self.essential_metabolites[i] = model.metabolites.get_by_id(m)

    def optimize_with_reaction(self, target, max_interventions=1, substrate=None,
                               biomass=None, design_method=OptKnock, max_results=100,
                               non_essential_metabolites=False, max_evaluations=20000, **design_kwargs):
        raise NotImplementedError

    def optimize_with_metabolites(self, target, max_interventions=1, substrate=None, biomass=None,
                                  max_results=100, non_essential_metabolites=False, max_evaluations=20000,
                                  **design_kwargs):
        raise NotImplementedError

    def essential_metabolites_reactions(self):
        essential_metabolites = self.model.essential_metabolites()
        reactions = set()
        for metabolite in essential_metabolites:
            reactions.update(metabolite.reactions)

        return reactions

    def _evaluate_designs(self, strain_designs, objective_function):
        evaluated_designs = DataFrame(columns=["design", "fitness"])

        for i, design in enumerate(strain_designs):
            with TimeMachine() as tm:
                design.apply(self.model, time_machine=tm)
                solution = fba(self.model)
                fitness = objective_function(self.model, solution, design.targets)
                evaluated_designs.loc[i] = [design, fitness]

        return evaluated_designs


class RandomMutagenesisDesign(GenericMARSIDesignMethod):
    """
    Apply only knockout like designs where total loss of functions are expected.



    """
    def optimize_with_reaction(self, target, max_interventions=1, substrate=None,
                               biomass=None, design_method="optgene", max_results=100,
                               non_essential_metabolites=False, max_evaluations=20000, **designer_kwargs):

        target_flux = self.model._reaction_for(target)

        exclude_reactions = []
        if non_essential_metabolites:
            exclude_reactions = self.essential_metabolites_reactions()

        if 'essential_reactions' in designer_kwargs:
            designer_kwargs['essential_reactions'] += exclude_reactions
        else:
            designer_kwargs['essential_reactions'] = exclude_reactions

        bpcy = biomass_product_coupled_yield(biomass, target_flux, substrate)

        if design_method == "optgene":
            designer = OptGene(self.model, **designer_kwargs)
            knockouts = designer.run(max_knockouts=max_interventions, biomass=biomass, substrate=substrate,
                                     target=target_flux, max_results=max_results, max_evaluations=max_evaluations,
                                     use_nullspace_simplification=False)
        elif design_method == "optknock":
            designer = OptKnock(self.model, **designer_kwargs)
            knockouts = designer.run(max_knockouts=max_interventions, biomass=biomass, substrate=substrate,
                                     target=target_flux, max_results=max_results)
        else:
            raise ValueError("'design_method' can be one of 'optgene' or 'optknock'")

        designs = self._evaluate_designs(iter(knockouts), bpcy)

        anti_metabolites_design = DataFrame()
        for i_, row in designs.iterrows():
            _anti_metabolites_design = replace_design(self.model, row.design, row.fitness, bpcy, fba, {},
                                                      ignore_metabolites=self.currency_metabolites,
                                                      essential_metabolites=self.essential_metabolites)
            anti_metabolites_design = anti_metabolites_design.append(_anti_metabolites_design, ignore_index=True)

        return anti_metabolites_design

    def optimize_with_metabolites(self, target, max_interventions=1, substrate=None, biomass=None,
                                  max_results=100, max_evaluations=20000, **design_kwargs):

        target_flux = self.model._reaction_for(target)

        designer = OptMet(model=self.model, essential_metabolites=CURRENCY_METABOLITES, **design_kwargs)
        knockouts = designer.run(max_knockouts=max_interventions, biomass=biomass, substrate=substrate,
                                 target=target_flux, max_results=max_results, max_evaluations=max_evaluations)

        return knockouts.data_frame


class ALEDesign(GenericMARSIDesignMethod):
    """
    Apply both knockout and flux modulation. The strains will be selected via growth rate proxy.


    """

    def optimize_with_reaction(self, target, max_interventions=1, substrate=None,
                               biomass=None, design_method="differential_fva", max_results=100,
                               non_essential_metabolites=False, **designer_kwargs):

        target_flux = self.model._reaction_for(target)
        bpcy = biomass_product_coupled_yield(biomass, target_flux, substrate)

        designer = DifferentialFVA(self.model, objective=target, variables=[biomass],
                                   points=max_results, **designer_kwargs)

        evaluated_designs = DataFrame(columns=["strain_designs", "fitness"])
        for i, design in enumerate(designer.run()):
            with TimeMachine() as tm:
                design.apply(self.model, time_machine=tm)
                solution = self.model.solve()
                fitness = bpcy(self.model, solution, design.targets)
                evaluated_designs.loc[i] = [design, fitness]

        anti_metabolites_design = DataFrame()
        for _, row in evaluated_designs.iterrows():
            _anti_metabolites_design = replace_design(self.model, row.design, row.fitness, bpcy, fba, {},
                                                      ignore_metabolites=self.currency_metabolites,
                                                      essential_metabolites=self.essential_metabolites)
            anti_metabolites_design = anti_metabolites_design.append(_anti_metabolites_design, ignore_index=True)

        return anti_metabolites_design

    def optimize_with_metabolites(self, target, max_interventions=1, substrate=None, biomass=None,
                                  max_results=100, max_evaluations=20000, **design_kwargs):
        raise NotImplementedError
