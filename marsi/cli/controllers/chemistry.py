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
from cement.core.controller import CementBaseController, expose

from marsi import config
from marsi.nearest_neighbors import search_closest_compounds
from marsi.chemistry.molecule import Molecule
from mongoengine import connect

from marsi.io import write_excel_file

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
            (['--atoms-diff', '-ad'], dict(help="The maximum number of the atoms difference for database query")),
            (['--bonds-diff', '-bd'], dict(help="The maximum number of the bonds difference for database query")),
            (['--rings-diff', '-rd'], dict(help="The maximum number of the rings difference for database query")),
            (['--output-file', '-o'], dict(help="Output file")),
            (['--output-format', '-f'], dict(help='Output format (default: csv)', default='csv'))
        ]

    @expose(hide=True)
    def default(self):
        print("Welcome to MARSI chemistry package")
        print("Here you can find the tools to find and sort analogs for metabolites")


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

        molecule = None

        bonds_weight = 0.5
        atoms_weight = 0.5
        bonds_diff = 3
        atoms_diff = 4
        rings_diff = 3

        if self.app.pargs.bonds_weight is not None:
            bonds_weight = float(self.app.pargs.bonds_weight)

        if self.app.pargs.atoms_weight is not None:
            atoms_weight = float(self.app.pargs.atoms_weight)

        if self.app.pargs.bonds_diff is not None:
            bonds_diff = float(self.app.pargs.bonds_diff)

        if self.app.pargs.atoms_diff is not None:
            atoms_diff = float(self.app.pargs.atoms_diff)

        if self.app.pargs.rings_diff is not None:
            rings_diff = float(self.app.pargs.rings_diff)

        if self.app.pargs.inchi is not None:
            molecule = Molecule.from_inchi(self.app.pargs.inchi)
        elif self.app.pargs.sdf is not None:
            molecule = Molecule.from_sdf(self.app.pargs.sdf)
        elif self.app.pargs.mol is not None:
            molecule = Molecule.from_mol(self.app.pargs.mol)
        else:
            print("Please provide one of the following inputs --inchi, --sdf or --mol")
            exit(1)

        fpformat = self.app.pargs.fingerprint_format

        results = search_closest_compounds(molecule=molecule, fpformat=fpformat, bonds_weight=bonds_weight,
                                           bonds_diff=bonds_diff, atoms_weight=atoms_weight, atoms_diff=atoms_diff,
                                           rings_diff=rings_diff)

        OUTPUT_WRITERS[self.app.pargs.output_format](results, output_file, None)
