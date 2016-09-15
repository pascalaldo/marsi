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

from pandas import melt, DataFrame

from bokeh.charts import Bar
from bokeh.io import show


def summary_plot(summary, dbs=['chebi', 'kegg', 'drugbank', 'pubchem']):
    summary = DataFrame(summary)
    summary['key'] = summary.index
    by_database = melt(summary, id_vars='key', value_vars=dbs).dropna(subset=['value'])
    bar = Bar(by_database, label='variable', values='key', color='variable', agg='count')
    show(bar)
