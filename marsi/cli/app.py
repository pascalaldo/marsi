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
from cement.core.foundation import CementApp

from marsi.cli.controllers import MarsiBaseController
from marsi.cli.controllers.chemistry import ChemistryController
from marsi.cli.controllers.database import DatabaseController
from marsi.cli.controllers.modeling import OptimizationController


class MarsiApp(CementApp):
    class Meta:
        label = 'marsi'
        base_controller = 'base'
        handlers = [
            MarsiBaseController,
            DatabaseController,
            OptimizationController,
            ChemistryController
        ]


def main():
    try:
        with MarsiApp() as app:
            app.run()
    except KeyboardInterrupt:
        print("Cancelled by user")
        exit(130)
