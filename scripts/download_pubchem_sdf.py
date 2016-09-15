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
import pubchempy as pbc
from IProgress.progressbar import ProgressBar
from IProgress.widgets import Bar, Percentage, ETA

from marsi.io.data import pubchem
from marsi.utils import data_dir

pubchem_ids = pubchem.compound_id.unique()
progress = ProgressBar(widgets=["Downloading: ", Bar(), Percentage(), " ", ETA()], maxval=len(pubchem_ids))

for pubchem_id in progress(pubchem_ids):
    try:
        pbc.download('sdf', os.path.join(data_dir, 'pubchem_sdf_files', '%i.sdf' % pubchem_id), int(pubchem_id))
    except IOError:
        # File already exists
        continue
