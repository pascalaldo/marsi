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


import cobra
import numpy

from networkx import all_shortest_paths, all_simple_paths
from networkx.classes.multidigraph import MultiDiGraph

from itertools import product


def tf_idf(model):
    """
    from Wikipedia (https://en.wikipedia.org/wiki/Tf–idf):
    "In information retrieval, tf–idf, short for term frequency–inverse document frequency, is a numerical statistic
    that is intended to reflect how important a word is to a document in a collection or corpus."

    This adapted version checks how important a metabolite is to a subsystem/pathway in model.

    Parameters
    ----------
    model: cobra.Model

    Returns
    -------
    np.matrix
        A matrix of n*m of n metabolites and m subsystems with the 'importance' of a metabolites in each subsystem.
    """
    assert isinstance(model, cobra.Model)
    subsystems = []
    for r in model.reactions:
        if r.subsystem not in subsystems:
            subsystems.append(r.subsystem)
    n_subsystems = len(subsystems)
    n_metabolites = len(model.metabolites)
    doc_count = numpy.zeros(n_metabolites)
    subsystem_count = numpy.zeros((n_metabolites, n_subsystems))

    for reaction in model.reactions:
        assert isinstance(reaction, cobra.Reaction)
        subsystem_index = subsystems.index(reaction.subsystem)
        for metabolite in reaction.metabolites:
            m_index = model.metabolites.index(metabolite)
            doc_count[m_index] += 1
            subsystem_count[m_index, subsystem_index] += 1

    occurrences = numpy.apply_along_axis(sum, 1, subsystem_count)
    tf = numpy.divide(subsystem_count.T, occurrences).T
    idf = n_subsystems/doc_count

    _tf_idf = numpy.apply_along_axis(numpy.multiply, 0, tf, idf)
    return _tf_idf


SKIP_DEFAULT = ['h_c', 'h_p', 'h_e', 'h2o_c', 'atp_c', 'pi_c', 'adp_c', 'h2o_p', 'ppi_c', 'nad_c', 'nadh_c', 'nadp_c',
                'pi_p', 'nadph_c', 'amp_c', 'co2_c', 'co2_p', 'co2_e', 'coa_c', 'ACP_c', 'nh4_c', 'nh4_p', 'nh4_e']


class GraphPathway(object):
    def __init__(self, model, nodes):
        assert isinstance(model, GraphBasedModel)
        self._model = model
        self._nodes = nodes

    def __eq__(self, other):
        if isinstance(other, GraphPathway):
            return self._nodes == other._nodes
        elif isinstance(other, list):
            return self._nodes == other
        else:
            return False

    @property
    def edges(self):
        edges = []
        for i in range(1, len(self._nodes)):
            data = self._model.graph.get_edge_data(self._nodes[i-1], self._nodes[i])
            edges.append(data)

        return edges

    def __repr__(self):
        parts = []
        for i, edge in enumerate(self.edges):
            via = "(" + ", ".join(["%s w:[%f]" % (k, v["weight"]) for k, v in edge.items()]) + ")"
            parts.append("%s -> %s ->" % (self._nodes[i].id, via))

        return " | ".join(parts) + " -> %s" % self._nodes[-1]


class GraphBasedModel(cobra.Model):
    """
    This graph-based model is a reaction to reaction network where the weight of the edges is determined
    by an adaptation of the TF-IDF normalization used in information retrieval.

    """

    def __init__(self, id_or_model=None, name=None, skip_metabolites=SKIP_DEFAULT):
        self.graph = MultiDiGraph()
        self.reactions = cobra.DictList()
        self.metabolites = cobra.DictList()
        self.subsystems = []
        self.name = name
        self.skip_metabolites = skip_metabolites
        self.tf_idf = None
        if isinstance(id_or_model, cobra.Model):
            self.id = id_or_model.id
            self.add_metabolites(id_or_model.metabolites)
            self.add_reactions(id_or_model.reactions)
        else:
            self.id = id_or_model

    def add_metabolites(self, metabolites):
        for m in metabolites:
            self.metabolites.append(m)

    def add_reactions(self, reactions):
        for reaction in reactions:
            assert isinstance(reaction, cobra.Reaction)

            if reaction.subsystem not in self.subsystems:
                self.subsystems.append(reaction.subsystem)

            for metabolite in reaction.metabolites:
                if metabolite not in self.metabolites:
                    self.add_metabolites([metabolite])

            self.reactions.append(reaction)
        self.tf_idf = tf_idf(self)
        self.populate_graph()

    def populate_graph(self):
        for reaction in self.reactions:
            subs_idx = self.subsystems.index(reaction.subsystem)
            _from = [m for m, c in reaction.metabolites.items() if c < 0]
            _to = [m for m, c in reaction.metabolites.items() if c > 0]

            for f, t in product(_from, _to):
                if t.id in self.skip_metabolites or f.id in self.skip_metabolites:
                    continue

                f_idx, t_idx = self.metabolites.index(f), self.metabolites.index(t)
                weight = (self.tf_idf.max() - self.tf_idf[t_idx][subs_idx]) +\
                         (self.tf_idf.max() - self.tf_idf[t_idx][subs_idx])

                if not self.graph.has_edge(f, t, reaction.id):
                    self.graph.add_edge(f, t, weight=weight, key=reaction.id)
                if reaction.reversibility:
                    if not self.graph.has_edge(t, f, reaction.id + "_rev"):
                        self.graph.add_edge(t, f, weight=weight, key=reaction.id + "_rev")

    def paths(self, metabolite_a, metabolite_b, radius=2):
        shortest_paths = list(all_shortest_paths(self.graph, metabolite_a, metabolite_b, weight="weight"))
        sizes = set([len(sp) for sp in shortest_paths])
        paths = [GraphPathway(self, sp) for sp in shortest_paths]
        simple_paths = all_simple_paths(self.graph, metabolite_a, metabolite_b, cutoff=max(sizes)+radius)
        for sp in simple_paths:
            if sp not in paths:
                paths.append(GraphPathway(self, sp))

        return paths
