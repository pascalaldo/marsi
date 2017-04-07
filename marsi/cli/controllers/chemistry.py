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
from numbers import Number

import os
import pybel
import six
from cement.core.controller import CementBaseController, expose
from pandas import DataFrame
from pymongo.errors import ServerSelectionTimeoutError

from marsi import config, utils
from marsi.nearest_neighbors import load_nearest_neighbors_model, build_feature_table
from marsi.chemistry import openbabel as ob, rdkit as rd
from mongoengine import connect

from marsi.io import write_excel_file
from marsi.io.mongodb import Metabolite, Database

from IProgress import ProgressBar, Percentage, ETA, Bar

from marsi.utils import pickle_large

OUTPUT_WRITERS = {
    'csv': lambda df, path, *args: df.to_csv(path),
    'excel': lambda df, path, *args: write_excel_file(df, path)

}


class ChemistryController(CementBaseController):
    """
    This is the Optimization Controller. It allows to run optimizations from the command line.

    """

    class Meta:
        label = 'chem'
        stacked_on = 'base'
        stacked_type = 'nested'
        description = "Metabolite database query and filter tools"
        arguments = [
            (['--inchi'], dict(help="The metabolite InChI to search")),
            (['--sdf'], dict(help="The metabolite SDF to search")),
            (['--mol'], dict(help="The metabolite MOL to search")),
            (['--fingerprint-format', '-fp'], dict(help="The fingerprint format", default='maccs', action='store')),
            (['--neighbors', '-k'], dict(help="Filter the first K hits")),
            (['--radius', '-r'], dict(help="Filter hits within R distance radius")),
            (['--atoms-weight', '-aw'], dict(help="The weight of the atoms for structural similarity")),
            (['--bonds-weight', '-bw'], dict(help="The weight of the bonds for structural similarity")),
            (['--output-file', '-o'], dict(help="Output file")),
            (['--output-format', '-f'], dict(help='Output format (default: csv)', default='csv'))
        ]

    @expose(hide=True)
    def default(self):
        print("Welcome to MARSI chemistry package")
        print("Here you can find the tools to find and sort analogs for metabolites")

    @expose(help="Build a nearest neighbors model with the specified fingerprints")
    def nearest_neighbors_model(self):
        connect(config.db_name)
        fpformat = self.app.pargs.fingerprint_format
        if fpformat not in pybel.fps:
            print("Invalid 'fingerprint-format' %s. Use of of %s" % (fpformat, ", ".join(pybel.fps)))

        solubility = "high"

        print("This process will take several hours, please confirm the following information:")
        print("+--------------------+---------------+")
        print("| Fingerprint format |        %s |" % fpformat.rjust(6))
        print("+--------------------+---------------+")
        print("|         Solubility |        %s |" % solubility.rjust(6))
        print("+--------------------+---------------+")

        value = input("Do you want to continue y/n [y]: ")

        if len(value) == 0:
            value = "y"

        if value != "y":
            print("Cancelled")
            exit(1)

        model_file = os.path.join(utils.data_dir, "fingerprints_default_%s_sol_%s.pickle" % (fpformat, solubility))
        print("File will be saved as %s" % model_file)

        try:
            indices, features, lengths = build_feature_table(Database.metabolites, fpformat=fpformat,
                                                             solubility=solubility, chunk_size=1e5)
            pickle_large((indices, features, lengths), model_file, progress=True)
            print("Completed")
        except ServerSelectionTimeoutError:
            print("Error: cannot connect to database. Make sure your mongodb is running and accessible.")
            exit(1)

    @expose(help="Find analogs for a metabolite")
    def find_analogs(self):
        """
        1. Make a fingerprint from the --inchi, --sdf or --mol.
        2. Query the database

        Returns
        -------

        """
        connect(config.db_name)
        output_file = None

        if self.app.pargs.output_file is not None:
            output_file = self.app.pargs.output_file
        else:
            print("--output-file argument is required")
            exit(1)

        output_file += ".%s" % self.app.pargs.output_format

        ob_molecule = None
        rd_molecule = None

        fingerprint = None
        fingerprint_db = None

        use_volume = False
        max_volume_percentage = 0.5

        bonds_weight = 0.5
        atoms_weight = 0.5

        if self.app.pargs.bonds_weight is not None:
            bonds_weight = float(self.app.pargs.bonds_weight)

        if self.app.pargs.atoms_weight is not None:
            atoms_weight = float(self.app.pargs.atoms_weight)

        if self.app.pargs.inchi is not None:
            ob_molecule = ob.inchi_to_molecule(self.app.pargs.inchi)
            rd_molecule = rd.inchi_to_molecule(self.app.pargs.inchi)
        elif self.app.pargs.sdf is not None:
            use_volume = True
            ob_molecule = ob.sdf_to_molecule(self.app.pargs.sdf)
            rd_molecule = rd.sdf_to_molecule(self.app.pargs.sdf)
        elif self.app.pargs.mol is not None:
            use_volume = True
            ob_molecule = ob.mol_to_molecule(self.app.pargs.mol)
            rd_molecule = rd.mol_to_molecule(self.app.pargs.mol)
        else:
            print("Please provide one of the following inputs --inchi, --sdf or --mol")
            exit(1)

        try:
            fpformat = self.app.pargs.fingerprint_format
            fingerprint = ob_molecule.calcfp(fpformat).fp
            fingerprint_db = load_nearest_neighbors_model(fpformat=fpformat)
        except ValueError as e:
            print(e)
            exit(1)

        if self.app.pargs.neighbors:
            neighbors = fingerprint_db.k_nearest_neighbors(fingerprint, int(self.app.pargs.neighbors))
        else:
            neighbors = fingerprint_db.radius_nearest_neighbors(fingerprint, float(self.app.pargs.radius or 0.75))

        print("Found %i neighbors" % len(neighbors))

        ob_descriptors = ob_molecule.calcdesc()
        descriptors = [key for key, value in six.iteritems(ob_descriptors) if isinstance(value, Number)]
        desc_keys = ["d_%s" % key for key in descriptors]

        ob_molecules = {}
        extra_descriptors = []

        ob_volume = None
        if use_volume:
            ob_volume = ob.monte_carlo_volume(ob_molecule)
            extra_descriptors.append('d_volume')

        results = DataFrame(columns=["tanimoto_similarity", "structural_similarity"] + desc_keys + extra_descriptors)

        progress = ProgressBar(maxval=len(neighbors), widgets=["Processing neighbors: ",
                                                               Bar(), Percentage(), "|", ETA()])
        progress.start()

        for inchi_key, tanimoto_distance in six.iteritems(neighbors):

            metabolite = Metabolite.get(inchi_key)

            ob_molecules[inchi_key] = _ob_molecule = metabolite.molecule('openbabel')
            _rd_molecule = metabolite.molecule('rdkit')
            extra_descriptors_values = []

            if use_volume:
                volume_diff = abs(ob_volume - ob.monte_carlo_volume(ob_molecules[inchi_key]))
                if volume_diff > max_volume_percentage * ob_volume:
                    progress.update(progress.currval + 1)
                    continue
                extra_descriptors_values.append(volume_diff)

            _ob_descriptors = _ob_molecule.calcdesc(descriptors)

            descriptors_changes = [abs(ob_descriptors[k] - _ob_descriptors[k]) for k in descriptors]
            structural_similarity = rd.structural_similarity(rd_molecule,
                                                             _rd_molecule,
                                                             bonds_weight=bonds_weight,
                                                             atoms_weight=atoms_weight)

            distances = [1 - tanimoto_distance, structural_similarity]

            results.loc[inchi_key] = distances + descriptors_changes + extra_descriptors_values
            progress.update(progress.currval + 1)
        progress.finish()

        OUTPUT_WRITERS[self.app.pargs.output_format](results, output_file, ob_molecules)
