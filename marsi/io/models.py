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

from IProgress.progressbar import ProgressBar
from IProgress.widgets import Percentage, Bar, ETA

from cameo import load_model as cameo_load_model
from marsi.utils import models_dir
from marsi.processing.models import find_inchi_for_bigg_metabolite
from cobra.io.json import save_json_model


def annotate_metabolite(metabolite):
    try:
        metabolite.annotation['inchi']
    except KeyError:
        try:
            metabolite.annotation['inchi'] = find_inchi_for_bigg_metabolite(metabolite.model, metabolite.id)
        except ValueError:
            pass


def annotate_model(model):
    pbar = ProgressBar(maxval=len(model.metabolites), widgets=["Annotating: ", Percentage(), Bar(), ETA()])
    [annotate_metabolite(m) for m in pbar(model.metabolites)]


def load_model(path_or_handle):
    if isinstance(path_or_handle, str) and os.path.isfile(os.path.join(models_dir, "%s.json" % path_or_handle)):
        return cameo_load_model(os.path.join(models_dir, "%s.json" % path_or_handle))

    model = cameo_load_model(path_or_handle)
    annotate_model(model)
    save_json_model(model, "%s.json" % model.id)
    return model
