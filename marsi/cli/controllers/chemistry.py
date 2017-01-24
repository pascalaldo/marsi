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
            (['--fp'], dict(help="The fingerprint format")),
            (['--k'], dict(help="Filter the first K hits")),
            (['--r'], dict(help="Filter hits within r distance radius"))
        ]

    @expose()
    def find_analogs(self):
        """
        1. Make a fingerprint from the --inchi, --sdf or --mol.
        2. Query the database

        Returns
        -------

        """
