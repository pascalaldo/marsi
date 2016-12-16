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
import traceback

import io
import os
import re

from cement.core.controller import CementBaseController, expose

from IProgress import ProgressBar, Bar, ETA

from cameo.io import load_model as load_model_for_optimization
from cameo.strain_design import OptGene

from marsi.design import GenericMARSIDesignMethod
from marsi.io.build_database import build_database
from marsi.io.parsers import parse_chebi_data, parse_kegg_brite, parse_pubchem
from marsi.io.retriaval import retrieve_drugbank_open_structures, retrieve_drugbank_open_vocabulary, \
    retrieve_bigg_reactions, retrieve_bigg_metabolites, retrieve_chebi_structures, retrieve_chebi_names, \
    retrieve_chebi_relation, retrieve_chebi_vertice, retrieve_kegg_brite, retrieve_zinc_properties, \
    retrieve_kegg_mol_files, retrieve_pubchem_mol_files
from marsi.utils import default_carbon_sources, unique, data_dir, log_dir, models_dir
from marsi.config import prj_dir

__all__ = ["MarsiBaseController"]


BIOMASS_RE = re.compile("biomass", re.IGNORECASE)


class MarsiBaseController(CementBaseController):
    """
    This is the application base controller.

    """

    class Meta:
        label = 'base'

    @expose(hide=True)
    def default(self):
        print("Welcome to MARSI. For more information type marsi --help")


class InitializationController(CementBaseController):
    download_steps = ["download_drug_bank",
                      "download_bigg",
                      "download_chebi",
                      "download_pubchem",
                      "download_zinc",
                      "download_kegg"]

    build_steps = ["build_chebi_data",
                   "build_pubchem_data",
                   "build_kegg_data",
                   "build_database"]

    class Meta:
        label = 'init'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Initialise MARSI (download data and build initial database)"
        arguments = [
            (['--drugbank-version'], dict(help="drugbank version (5.0.3)"))
        ]

    @expose(help="Init MARSI", hide=True)
    def default(self):
        """
        1. Download necessary files.
        2. Build the data files.
        3. Build the database.

        Returns
        -------

        """
        init_file_path = os.path.join(prj_dir, ".init")
        in_progress = os.path.exists(init_file_path)

        if not os.path.exists(data_dir):
            os.mkdir(data_dir)
            os.mkdir(log_dir)
            os.mkdir(models_dir)

        i, init_file = self._check_state("download", init_file_path, in_progress)

        self._run_sequence(self.download_steps, i, "Downloading files: ", "Download complete", init_file)

        i, init_file = self._check_state("build", init_file_path, in_progress)

        self._run_sequence(self.download_steps, i, "Building database: ", "Build complete", init_file)
        init_file.close()

    def _run_sequence(self, steps, i, pbar_label, complete_message, state_file):
        pbar = ProgressBar(maxval=len(self.build_steps), widgets=[pbar_label, Bar(), ETA()])
        pbar.start()
        pbar.update(i-1)

        for i in range(i, len(self.download_steps) - 1):
            method = getattr(self, self.download_steps[i])
            try:
                method()
                state_file.write(self.download_steps[i])
                pbar.update(i)
            except Exception as e:
                state_file.close()
                traceback.print_exc()
                print("Something went wrong while running %s: %s" % (self.download_steps[i], str(e)))
                exit(1)

        pbar.finish()
        print(complete_message)

    def _check_state(self, label, state_file_path, in_progress):
        i = 0
        if in_progress:
            init_file = open(state_file_path, "a+")
            for line in init_file:
                if line.startswith(label):
                    assert self.download_steps[i] == line
                    i += 1
                else:
                    break
        else:
            init_file = open(state_file_path, "w")

        return i, init_file

    def download_drug_bank(self):
        retrieve_drugbank_open_structures(self.app.pargs.drugbank_version or "5.0.3")
        retrieve_drugbank_open_vocabulary(self.app.pargs.drugbank_version or "5.0.3")

    @staticmethod
    def download_bigg():
        retrieve_bigg_reactions()
        retrieve_bigg_metabolites()

    @staticmethod
    def download_chebi():
        retrieve_chebi_structures()
        retrieve_chebi_names()
        retrieve_chebi_relation()
        retrieve_chebi_vertice()

    def download_pubchem(self):
        pass

    @staticmethod
    def download_zinc():
        retrieve_zinc_properties()

    @staticmethod
    def download_kegg():
        retrieve_kegg_brite()

    @staticmethod
    def build_chebi_data():
        chebi_names_file = os.path.join(data_dir, "chebi_names_3star.txt")
        chebi_vertice_file = os.path.join(data_dir, "chebi_vertice_3star.tsv")
        chebi_relation_file = os.path.join(data_dir, "chebi_relation_3star.tsv")
        chebi_data = parse_chebi_data(chebi_names_file, chebi_vertice_file, chebi_relation_file)
        chebi_data.to_tsv(os.path.join(data_dir, "chebi_analogues_filtered.csv"))

    @staticmethod
    def build_pubchem_data():
        pubchem = parse_pubchem(os.path.join(data_dir, "pubchem_compound_analogs_antimetabolites.txt"))
        pubchem.to_csv(os.path.join(data_dir, "pubchem_data.csv"))
        retrieve_pubchem_mol_files(pubchem.compound_ids.unique())

    @staticmethod
    def build_kegg_data():
        kegg = parse_kegg_brite(os.path.join(data_dir, "kegg_brite_08310.keg"))
        kegg.to_csv(os.path.join(data_dir, "kegg_data.csv"))
        retrieve_kegg_mol_files(kegg, dest=data_dir)

    @staticmethod
    def build_database():
        from marsi.io import data
        build_database(data, data_dir)


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
            (['--model'], dict(help="path or identifier of the model to be used")),
            (['--carbon-source'], dict(help="(optional) the carbon source exchange reaction. "
                                            "It will be auto-detected if not defined")),
            (['--target'], dict(help="The exchange reaction of the target phenotype or the metabolite to accumulate")),
            (['--approach'], dict(help="classic: use reaction search; metabolites: search for "
                                       "metabolite targets directly. (default: classic)")),
            (['--biomass'], dict(help="(optional) the biomass reaction (needed if there is none or many reactions "
                                      "containing biomass on their ids or names")),
            (['--max-interventions', '-m'], dict(help="(optional) maximum number of interventions (default; 2)")),
            (['--output-file', '-o'], dict(help="output file")),
            (['--transporters'], dict(help="Include transporters knockout"))
        ]

    @expose(hide=True)
    def default(self):
        """
        1. Load model specified by --model.
        2. Detect carbon source if not defined, otherwise change the carbon source.
        3. Run the optimization.
        4. Return the results.

        Returns
        -------

        """
        target = None
        model = None
        carbon_source = None
        max_interventions = 2
        output_file = None

        if self.app.pargs.output_file is None:
            print("--output-file is required")
            exit(1)
        else:
            output_file = self.app.pargs.output_file

        if self.app.pargs.max_interventions is not None:
            max_interventions = int(self.app.pargs.max_interventions)

        if self.app.pargs.model is None:
            print("--model is required")
            exit(1)
        else:
            try:
                model = load_model_for_optimization(self.app.pargs.model)
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
                print("Please chose:")
                print("\n".join(str(s.id) for s in biomass))
                exit(1)
            else:
                biomass = biomass[0]

        designer = GenericMARSIDesignMethod(model=model)

        if self.app.pargs.transporters is None:
            transporters = [r.id for r in model.reactions if len(set(m.compartment for m in r.metabolites)) > 1]
            exclude_reactions = self.exclude_reactions
            exclude_reactions += transporters
        else:
            exclude_reactions = self.exclude_reactions

        if self.app.pargs.approach == 'classic':
            result = designer.optimize_with_reaction(target=target, substrate=carbon_source, biomass=biomass,
                                                     design_method=OptGene, manipulation_type="reactions",
                                                     max_interventions=max_interventions,
                                                     essential_reactions=exclude_reactions)
        else:
            result = designer.optimize_with_metabolites(target=target, substrate=carbon_source, biomass=biomass,
                                                        max_interventions=max_interventions)

        result.to_csv(output_file)

        print("Finished")
