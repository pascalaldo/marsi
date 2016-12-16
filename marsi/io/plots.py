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

from itertools import combinations

from pandas import melt, DataFrame

from bokeh.charts import Bar
from bokeh.io import show


def summary_plot(summary, dbs=['chebi', 'kegg', 'drugbank', 'pubchem']):
    data = DataFrame(summary)
    data['key'] = data.index
    by_database = melt(data, id_vars='key', value_vars=dbs).dropna(subset=['value'])
    bar = Bar(by_database, label='variable', values='key', color='variable', agg='count')
    show(bar)


def venn_plot(summary, dbs=['chebi', 'kegg', 'drugbank', 'pubchem']):
    sets = {db: set(summary[~summary[db].isnull()].index) for db in dbs}
    for db1, db2 in combinations(dbs, 2):
        l_int = len(sets[db1].intersection(sets[db2]))
        print("%s AND %s: %i" % (db1, db2, l_int))

    for db1, db2, db3 in combinations(dbs, 3):
        l_int = len(sets[db1].intersection(sets[db2]).intersection(db3))
        print("%s AND %s AND %s: %i" % (db1, db2, db3, l_int))

    for db1, db2, db3, db4 in combinations(dbs, 4):
        l_int = len(sets[db1].intersection(sets[db2]).intersection(db3).intersection(db4))
        print("%s AND %s AND %s AND %s: %i " % (db1, db2, db3, db4, l_int))
