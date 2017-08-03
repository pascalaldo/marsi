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
import pytest
from cameo.flux_analysis.analysis import n_carbon

from marsi.utils import frange, default_carbon_sources, unique


def test_frange():
    range_with_start_and_stop = frange(1.1, 2.1, steps=10)
    value = 1.1
    assert next(range_with_start_and_stop) == pytest.approx(value, 1e-6)
    has_next = True
    while has_next:
        try:
            next_val = next(range_with_start_and_stop)
            value += 0.1
            assert next_val == pytest.approx(value, 1e-6)
        except StopIteration:
            has_next = False

    assert value == pytest.approx(2, 1e-4)

    range_without_stop = frange(4, steps=20)

    value = 0
    assert next(range_without_stop) == pytest.approx(value, 1e-6)
    has_next = True
    while has_next:
        try:
            next_val = next(range_without_stop)
            value += 0.2
            assert next_val == pytest.approx(value, 1e-6)
        except StopIteration:
            has_next = False

    assert value == pytest.approx(3.8, 1e-4)


def test_default_carbon_sources(model):
    carbon_sources = default_carbon_sources(model)
    assert len(carbon_sources) == 1
    assert n_carbon(carbon_sources[0]) > 0
    assert carbon_sources[0] == model.reactions.EX_glc__D_e


def test_unique():
    list_a = [1, 1, 3, 5, 5, 7, 9]
    unique(list_a)

    assert list_a == [1, 3, 5, 7, 9]
