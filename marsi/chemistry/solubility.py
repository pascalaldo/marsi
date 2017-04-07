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

import numpy as np
from scipy.optimize import curve_fit
from marsi.processing.chemistry.openbabel import smiles_to_molecule


def rss(measured, predicted):
    assert isinstance(measured, np.ndarray)
    assert isinstance(predicted, np.ndarray)

    residuals = (measured - predicted)**2
    return sum(residuals)


def fit_model(measured, smiles):
    molecules = [smiles_to_molecule(_smiles) for _smiles in smiles]
    return curve_fit(calc_solubility, molecules, measured)


def calc_solubility(molecules, i, log_p, mw, hbd, hba, mp, mr, abonds, rb, ap):
    y = [0 for _ in molecules]

    for j, molecule in enumerate(molecules):
        desc = molecule.calcdesc()

        rotatable_bounds = molecule.OBMol.NumRotors()
        aromatic_partition = len([a for a in molecule.atoms if a.OBAtom.IsAromatic()]) / desc['atoms']

        y[j] = i + log_p * desc['logP'] + mw * desc['MW'] + hbd * desc['HBD'] + hba * desc['HBA2'] + \
               mp * desc['MP'] + mr * desc['MR'] + abonds * desc['abonds'] + rb * rotatable_bounds + \
               ap * aromatic_partition

    return np.array(y)
