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

from pandas import read_csv, DataFrame

from marsi.utils import internal_data_dir

try:
    bigg_metabolites = read_csv(os.path.join(internal_data_dir, "bigg_models_metabolites.txt"), sep="\t")
    bigg_metabolites.database_links = bigg_metabolites.database_links.apply(eval)
    bigg_metabolites.model_list = bigg_metabolites.model_list.apply(str.split, args=(", ",))
except IOError:
    bigg_metabolites = DataFrame()
