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
"""
Quantitative Structure-Activity Relationship.

This module contains basic Machine Learning stuff to analyse the IC50 data:
    1. Feature extraction;
    2. Model training;
    3. Feature analysis;
"""
import numpy
from pandas import Series, DataFrame

from marsi.processing.chemistry.openbabel import get_spectrophore_data


def _build_feature_row(metabolite, *functions, **params):
    return Series(name=metabolite.inchi_key,
                  data=[f(metabolite, **params) for f in functions],
                  index=[f.__name__ for f in functions])


def build_custom_feature_table(metabolites, *functions, **params):
    result = DataFrame()

    for metabolite in metabolites:
        result = result.append(_build_feature_row(metabolite, *functions, **params))

    return result


def build_target_table(ic50s):
    result = DataFrame(columns=["inchi_key", "target_name", "target_organism", "ic50"])
    for i, ic50 in enumerate(ic50s):
        result.loc[i] = [ic50.metabolite.inchi_key, ic50.target_name, ic50.target_organism, ic50.value]

    return result


def build_spectrophore_feature_table(metabolites):
    result = DataFrame(columns=["spectrophore_feature_%i" % (i + 1) for i in range(48)])
    for metabolite in metabolites:
        mol = metabolite.molecule(library='openbabel')
        mol.make3D()
        result.loc[metabolite.inchi_key] = get_spectrophore_data(mol)

    return result


def build_maccs_feature_table(metabolites):
    result = DataFrame(columns=["maccs%i" % (i + 1) for i in range(165)])
    for metabolite in metabolites:
        hashed = metabolite.fingerprint('maccs', hash=True, bits=165)
        result.loc[metabolite.inchi_key] = numpy.array(hashed, dtype=int)

    return result


def build_table(ic50s, use_maccs=True, use_spectrophore=True, *functions, **kwargs):
    target_table = build_target_table(ic50s)
    metabolites = set(ic50.metabolite for ic50 in ic50s)
    if use_maccs:
        print("Building maccs feature table")
        maccs_table = build_maccs_feature_table(metabolites)
        target_table = target_table.merge(maccs_table, right_index=True, left_on='inchi_key')
    if use_spectrophore:
        print("Building spectrophore feature table")
        spectrophore_table = build_spectrophore_feature_table(metabolites)
        target_table = target_table.merge(spectrophore_table, right_index=True, left_on='inchi_key')
    if len(functions) > 0:
        print("Building custom feature table (%in functions)" % len(functions))
        custom_table = build_custom_feature_table(metabolites, *functions, **kwargs)
        target_table = target_table.merge(custom_table, right_index=True, left_on='inchi_key')

    return target_table
