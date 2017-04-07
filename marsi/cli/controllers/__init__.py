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
import os
from cement.core.controller import CementBaseController, expose
from mongoengine import connect

from marsi.config import prj_dir, db_name
from marsi.utils import data_dir, log_dir, models_dir


class MarsiBaseController(CementBaseController):
    """
    This is the application base controller.

    """

    class Meta:
        label = 'base'

    @expose(hide=True)
    def default(self):
        """
        1. Download necessary files.
        2. Build the data files.
        3. Build the database.

        Returns
        -------

        """
        if not os.path.isdir(data_dir):
            os.mkdir(data_dir)
        if not os.path.isdir(log_dir):
            os.mkdir(log_dir)
        if not os.path.isdir(models_dir):
            os.mkdir(models_dir)
        print("\n################################################################")
        print("#                                                              #")
        print("#                       Welcome to MARSI!                      #")
        print("#                                                              #")
        print("################################################################")

        print("Your current project directory is %s" % prj_dir)
        print("The directory tree is:")
        print("%s:" % prj_dir)
        print("-- Data files: %s" % data_dir.split(os.sep)[-1])
        print("-- Log files:  %s" % log_dir.split(os.sep)[-1])
        print("-- Models dir: %s" % models_dir.split(os.sep)[-1])
        print("Your current database is %s" % db_name)
        print("\n")
        print("If you want to change this settings you can create a 'setup.cfg' as follows:\n")
        print("######### Begin MARSI configuration #######")
        print("[marsi]                                    ")
        print("# Working directory                        ")
        print("prj_dir=<directory of your choice>         ")
        print("# Name of the database                     ")
        print("db_name='marsi-db'                         ")
        print("######### End of MARSI configuration ######\n")

    @expose(help="Show db status")
    def db_status(self):
        from marsi.io.mongodb import Database

        connect(db_name)
        print("Database status (%s):" % db_name)
        print("|----------------------------------------------|--------------|")
        print("| Collection                                   | Total        |")
        print("|----------------------------------------------|--------------|")
        print("| Chemical compounds                           | %s |" % str(len(Database.metabolites)).ljust(12))
        print("|----------------------------------------------|--------------|")
        print("| Half inhibition concentrations (IC50)        | %s |" % str(len(Database.ic50s)).ljust(12))
        print("|----------------------------------------------|--------------|")
        print("| Half maximal effective concentrations (EC50) | %s |" % str(len(Database.ic50s)).ljust(12))
        print("|----------------------------------------------|--------------|")
