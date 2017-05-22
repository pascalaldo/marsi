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
from cameo import pfba, load_model
from cameo.strain_design.heuristic.evolutionary.objective_functions import biomass_product_coupled_yield
from cement.core.controller import CementBaseController, expose
from pandas import DataFrame

from marsi.utils import BIOMASS_RE
from marsi.cobra.strain_design import RandomMutagenesisDesign, ALEDesign
from marsi.cobra.utils import CURRENCY_METABOLITES
from marsi.cobra.strain_design.post_processing import replace_design
from marsi.utils import default_carbon_sources, unique


class OptimizationController(CementBaseController):
    """
    This is the Optimization Controller. It allows to run optimizations from the command line.

    """
    exclude_reactions = ["ATPM"]

    class Meta:
        label = 'optimize'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Optimize a host phenotype using Metabolites as targets"
        arguments = [
            (['--model', '-m'], dict(help="path or identifier of the model to be used")),
            (['--carbon-source'], dict(help="(optional) the carbon source exchange reaction. "
                                            "It will be auto-detected if not defined")),
            (['--target', '-t'], dict(help="The exchange reaction of the target phenotype or "
                                           "the metabolite to accumulate")),
            (['--approach', '-a'], dict(help="classic: use reaction search; metabolites: search for "
                                             "metabolite targets directly. (default: classic)")),
            (['--biomass'], dict(help="(optional) the biomass reaction (needed if there is none or many reactions "
                                      "containing biomass on their ids or names")),
            (['--max-interventions', '-mi'], dict(help="(optional) maximum number of interventions (default: 2)")),
            (['--max-evaluations', '-me'], dict(help='(optional) maximum number of evaluations (default: 20000)')),
            (['--output-file', '-o'], dict(help="output file")),
            (['--input-file', '-i'], dict(help="input file")),
            (['--transporters'], dict(help="Include transporters knockout"))
        ]

    @expose(hide=True)
    def default(self):
        pass

    @expose("Run optimization assuming Chemical Mutagenesis scenario")
    def mutagenesis(self):
        """
        1. Load model specified by --model.
        2. Detect carbon source if not defined, otherwise change the carbon source.
        3. Detect biomass if not defined
        4. Run the optimization.
        5. Return the results.
        """

        model, target, carbon_source, biomass, output_file = self._parse_optimization_args()

        designer = RandomMutagenesisDesign(model=model)

        max_interventions = 2
        max_evaluations = 20000

        if self.app.pargs.max_interventions is not None:
            max_interventions = int(self.app.pargs.max_interventions)

        if self.app.pargs.max_evaluations is not None:
            max_evaluations = int(self.app.pargs.max_evaluations)

        if self.app.pargs.transporters is not None:
            transporters = [r.id for r in model.reactions if len(set(m.compartment for m in r.metabolites)) > 1]
            exclude_reactions = self.exclude_reactions
            exclude_reactions += transporters
        else:
            exclude_reactions = self.exclude_reactions

        print("---------------------------------------------------")
        print("Chemical Mutagenesis Optimization setup:")
        print("- Search method: %s" % self.app.pargs.approach)
        print("- Target: %s" % target.id)
        print("- Biomass: %s" % biomass.id)
        print("- Carbon source: %s" % carbon_source.id)
        print("- Max evaluations: %i" % max_evaluations)
        print("- Max interventions: %i" % max_interventions)

        if self.app.pargs.approach == 'classic':
            result = designer.optimize_with_reaction(target=target, substrate=carbon_source, biomass=biomass,
                                                     design_method="optgene", manipulation_type="reactions",
                                                     max_interventions=max_interventions,
                                                     essential_reactions=exclude_reactions,
                                                     max_evaluations=max_evaluations)
        else:
            result = designer.optimize_with_metabolites(target=target, substrate=carbon_source, biomass=biomass,
                                                        max_interventions=max_interventions,
                                                        max_evaluations=max_evaluations)

        result.to_csv(output_file)

        print("Finished")

    @expose(help="Run optimization assuming Adaptive Laboratory Evolution scenario")
    def ale(self):
        """
                1. Load model specified by --model.
                2. Detect carbon source if not defined, otherwise change the carbon source.
                3. Detect biomass if not defined
                4. Run the optimization.
                5. Return the results.
                """

        model, target, carbon_source, biomass, output_file = self._parse_optimization_args()

        designer = ALEDesign(model=model)

        max_interventions = 2
        max_evaluations = 20000

        if self.app.pargs.max_interventions is not None:
            max_interventions = int(self.app.pargs.max_interventions)

        if self.app.pargs.max_evaluations is not None:
            max_evaluations = int(self.app.pargs.max_evaluations)

        if self.app.pargs.transporters is not None:
            transporters = [r.id for r in model.reactions if len(set(m.compartment for m in r.metabolites)) > 1]
            exclude_reactions = self.exclude_reactions
            exclude_reactions += transporters
        else:
            exclude_reactions = self.exclude_reactions

        print("---------------------------------------------------")
        print("ALE Optimization setup:")
        print("- Search method: %s" % self.app.pargs.approach)
        print("- Target: %s" % target.id)
        print("- Biomass: %s" % biomass.id)
        print("- Carbon source: %s" % carbon_source.id)
        print("- Max evaluations: %i" % max_evaluations)
        print("- Max interventions: %i" % max_interventions)

        if self.app.pargs.approach == 'classic':
            result = designer.optimize_with_reaction(target=target, substrate=carbon_source, biomass=biomass,
                                                     design_method="differential_fva", manipulation_type="reactions",
                                                     max_interventions=max_interventions,
                                                     essential_reactions=exclude_reactions)
        else:
            result = designer.optimize_with_metabolites(target=target, substrate=carbon_source, biomass=biomass,
                                                        max_interventions=max_interventions,
                                                        max_evaluations=max_evaluations)

        result.to_csv(output_file)

        print("Finished")

    @expose(help="Convert output from OptGene or OptKnock to Anti-metabolite approach")
    def convert(self):
        """
        1. Load model specified by --model.
        2. Detect carbon source if not defined, otherwise change the carbon source.
        3. Detect biomass if not defined
        4. Convert the results.

        Returns
        -------

        """
        input_file = None

        if self.app.pargs.input_file is None:
            print("--input-file is required")
            exit(1)
        else:
            input_file = self.app.pargs.input_file

        knockouts = DataFrame.from_csv(input_file)
        model, target, carbon_source, biomass, output_file = self._parse_optimization_args()

        objective_function = biomass_product_coupled_yield(biomass, target, carbon_source)

        knockouts['designs'] = knockouts.reactions.apply(lambda s: s)

        result = replace_design(model, knockouts, objective_function, pfba, {}, CURRENCY_METABOLITES)
        result.to_csv(output_file)

        print("Finished")

    def _parse_optimization_args(self):
        target = None
        model = None
        carbon_source = None
        output_file = None

        if self.app.pargs.output_file is None:
            print("--output-file is required")
            exit(1)
        else:
            output_file = self.app.pargs.output_file

        if self.app.pargs.model is None:
            print("--model is required")
            exit(1)
        else:
            try:
                model = load_model(self.app.pargs.model)
            except Exception:
                print("Invalid model %s" % self.app.pargs.model)
                exit(1)

        if self.app.pargs.target is None:
            print("--target is required")
            exit(1)
        else:
            try:
                target = model._reaction_for(self.app.pargs.target)
            except KeyError:
                print("Invalid target %s" % self.app.pargs.target)
                exit(1)

        if self.app.pargs.approach is None:
            self.app.pargs.approach = 'classic'

        model_carbon_source = default_carbon_sources(model)[0]
        if self.app.pargs.carbon_source is not None:
            try:
                carbon_source = model._reaction_for(self.app.pargs.carbon_source)
            except KeyError:
                print("Invalid carbon source %s" % self.app.pargs.carbon_source)
                exit(1)

            carbon_source.lower_bound = model_carbon_source.lower_bound
            model_carbon_source.lower_bound = 0
        else:
            carbon_source = model_carbon_source

        if self.app.pargs.biomass is not None:
            biomass = model.reactions.get_by_id(self.app.pargs.biomass)
        else:
            def is_biomass(name):
                if len(BIOMASS_RE.findall(name)) > 0:
                    return True
                else:
                    return False

            biomass = list(model.reactions.query(is_biomass, "id"))
            biomass += list(model.reactions.query(is_biomass, "name"))
            unique(biomass)

            if len(biomass) == 0:
                print("Cannot find biomass reaction")
                exit(1)
            elif len(biomass) > 1:
                print("Multiple biomass reactions found!")
                print("\n".join("%i: %s" % (i + 1, b.id) for i, b in enumerate(biomass)))
                choice = input("Please chose: ")
                try:
                    choice = int(choice)
                    biomass = biomass[choice - 1]
                except (IndexError, TypeError):
                    print("Invalid choice %i" % choice)
                    exit(1)
            else:
                biomass = biomass[0]

        return model, target, carbon_source, biomass, output_file
