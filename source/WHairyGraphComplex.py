""" !!!! DEPRECATED !!!!!!
Rather than this class, always use WRHairyGraphComplex (...it is quasi-isomorphic but smaller.)


Graph complexes based on simple graphs with numbered hairs and hairs of two decorations.
Implemented Differentials: Contract edges.

Graphs are realized as simple graphs with 2 extra vertices for epsilon and omega, (index self.n_vertices and self.n_vertices+1).
WARNING: If there is a tadpole the corresponding loop is not part of the graph--one can determine the presence of the tadpole
from the overall one too small loop number.
TODO: Take care that this does not produce problems
"""


__all__ = ['graph_type', 'sub_types', 'HairyGraphVS', 'HairyGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
           'EdgeToOneHairGO', 'EdgeToOneHairD', 'HairyGC']

import itertools
from sage.all import *
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import NautyInterface
import OrdinaryGraphComplex
import StoreLoad
import Parameters

graph_type = "whairy"

# Option to include zero hairs in the hairy graph complexes.
zero_hairs = False


def dict_to_list(d, n):
    return [(d[j] if j in d else j) for j in range(n)]


# ------- Graph Vector Space --------
class WHairyGraphVS(GraphVectorSpace.GraphVectorSpace):
    """Hairy graph vector space.

    Sub vector space with specified number of vertices, loops, hairs, even or odd edges, even or odd hair vertices
    and at least trivalent vertices. No multiple edges and not more than one hair is attached to a vertex. One hair is
    composed of a hair vertex and an edge connecting it to a vertex. The parity of the hair refers to the parity of the
    hair vertex alone.

    Attributes:
        - n_vertices (int): Number of internal vertices.
        - n_loops (int): Number of loops.
        - n_hairs (int): Number of hairs.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs (bool): Parity of the hair vertices. True for even hairs, False for odd hairs.
        - n_edges (int): Number of edges.
        - sub_type (str): Sub type of graphs.
        - ogvs (OrdinaryGraphComplex.OrdinaryGVS): Ordinary graph vector space without hairs.

    """

    def __init__(self, n_vertices, n_loops, n_hairs, n_ws):
        """Initialize the hairy graph vector space.

        :param n_vertices: Number of internal vertices. TODO: sure it is only internal?
        :type n_vertices: int
        :param n_loops: the genus of the graph.
        :type n_loops: int
        :param n_hairs: Number of hairs. They are distinguishable, numbered 1,...,n
        :type n_hairs: int
        :param n_ws: Number of w decorated hairs.
        :type even_edges: int
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs = n_hairs
        self.n_ws = n_ws
        # we count only the internal edges and omega and eps edges, but not the hair edges
        self.n_edges = self.n_loops + self.n_vertices
        self.sub_type = "w"
        super().__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(
            self.n_vertices + self.n_hairs+2, self.n_loops, False)

    def get_type(self):
        return 'wgraphs'

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs and self.n_ws == other.n_ws

    def __hash__(self):
        return hash("wgra%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wgra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "wgra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('hairs', self.n_hairs), ('ws', self.n_ws)])

    def get_partition(self):
        # All internal vertices are in color 1, the single eps vertex in color 2, the w vertex in color 3
        # and the hair vertices are in colors 4,...,n+3.
        return [list(range(self.n_vertices))] + [[j] for j in range(self.n_vertices, self.n_vertices + self.n_hairs+2)]

    def plot_graph(self, G):
        GG = Graph(G, loops=True)
        # add proper tadpole if needed
        if GG.size() < self.n_edges + self.n_hairs:
            GG.add_edge(self.n_vertices, self.n_vertices)

        return GG.plot(partition=self.get_partition(), vertex_labels=True)

    def is_valid(self):
        # At least trivalent internal vertices.
        l = (3 * self.n_vertices + self.n_ws <=
             2 * self.n_edges + self.n_hairs)
        # Nonnegative number of vertices, non negative number of loops, non-negative or positive number of hairs,
        # and at least one omega hair.
        l = l and self.n_vertices >= 0 and self.n_loops >= 0 and self.n_hairs >= 0 and self.n_ws >= 1
        # At most a full graph.
        l = l and self.n_edges <= (
            self.n_vertices+2) * (self.n_vertices + 1) / 2
        return l

    def get_work_estimate(self):
        # TODO
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_hairy_graphs(self, nvertices, nloops, nhairs, include_novertgraph=false):
        """ Produces all connected hairy graphs with nhairs hairs, that are the last vertices in the ordering.
        Graphs can have multiple hairs, but not tadpoles or multiple edges.
        :param include_novertgraph: Whether to include the graph with one edge and no vertices as a two-hair graph
        :type include_novertgraph: bool
        """
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2.
        # Degree 1 vertices are hairs, degree 2 vertices are edges and are removed later.
        nedges = nloops + nvertices - 1  # number of internal edges
        n_vertices_1 = nvertices
        n_vertices_2 = nhairs + nedges
        n_edges_bip = nhairs + 2 * nedges
        deg_range_1 = (3, n_edges_bip + 1)
        deg_range_2 = (1, 2)
        bipartite_graphs = NautyInterface.list_bipartite_graphs2(
            n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip)
        unordered = [self._bip_to_ordinary(
            G, nvertices, nedges, nhairs) for G in bipartite_graphs]
        # Produce all permutations of the hairs
        #all_perm = [ range(0,nvertices) + p for p in Permutations(range(nvertices, nvertices+nhairs)) ]
        # return [G.relabel(p, inplace=False) for p in all_perm ]
        if include_novertgraph and nvertices == 0 and nhairs == 2 and nloops == 0:
            unordered.append(Graph([(0, 1)]))
        return unordered

    def get_generating_graphs(self):
        # Idea: produce all hairy graphs with either two additional vertices or
        # one additional vertex , or an additional vertex and a hair.
        # hair, the second color being either of degree 1 or 2.
        # Degree 1 vertices are hairs, degree 2 vertices are edges and are removed later.
        # Depending on the total number of extra hairs (omega or eps) one makes a color 1 vertex (if >=3),
        # or a bivalent or univalent color 2 vertex into the extra hairs
        # No multiple hairs and edges.
        if not self.is_valid():
            return []

        if self.n_ws >= 2:
            # graphs w/o omega tadpole, and valency of the omega and eps vertex at least two
            hairy = self.get_hairy_graphs(
                self.n_vertices+2, self.n_loops, self.n_hairs)
            ret1 = [GG for G in hairy for GG in self._hairy_to_w_edge(G)]
            # graphs with omega-tadpole, and valency of the extra vertices at least 2, not counting the tadpole
            hairy = self.get_hairy_graphs(
                self.n_vertices+2, self.n_loops-1, self.n_hairs)
            ret2 = [GG for G in hairy for GG in self._hairy_to_w_edge(
                G, addtadpole=True)]
            # graphs w/o omega tadpole, and valency of the eps vertex equals 1
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops, self.n_hairs)
            ret3 = [GG for G in hairy for GG in self._hairy_to_w_edge_1eps(G)]
            # graphs with omega-tadpole, and valency of the eps vertex 1, not counting the tadpole
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops-1, self.n_hairs)
            ret4 = [GG for G in hairy for GG in self._hairy_to_w_edge_1eps(
                G, addtadpole=True)]
            # graphs w/o omega tadpole, and valency of the eps vertex equals 0
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops, self.n_hairs+1)
            ret5 = [GG for G in hairy for GG in self._hairy_to_w_hair_2w(G)]
            # graphs with omega-tadpole, and valency of the eps vertex 0, not counting the tadpole
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops-1, self.n_hairs+1)
            ret6 = [GG for G in hairy for GG in self._hairy_to_w_hair_2w(
                G, addtadpole=True)]

            # graphs w/o omega tadpole, and valency of the eps vertex equals 1
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops-1, self.n_hairs+2)
            ret7 = [
                GG for G in hairy for GG in self._hairy_to_w_hairhair_2w(G)]
            # graphs with omega-tadpole, and valency of the eps vertex 1, not counting the tadpole
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops-2, self.n_hairs+2)
            ret8 = [GG for G in hairy for GG in self._hairy_to_w_hairhair_2w(
                G, addtadpole=True)]
            ret9 = []
            ret10 = []
        else:
            # graphs w/o omega tadpole, and valency of the eps vertex at least two
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops, self.n_hairs)
            ret1 = [GG for G in hairy for GG in self._hairy_to_w_edge_1w(G)]
            # graphs with omega-tadpole, and valency of the eps vertex at least 2, not counting the tadpole
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops-1, self.n_hairs)
            ret2 = [GG for G in hairy for GG in self._hairy_to_w_edge_1w(
                G, addtadpole=True)]
            # graphs w/o omega tadpole, and valency of the eps vertex equals 1
            hairy = self.get_hairy_graphs(
                self.n_vertices, self.n_loops, self.n_hairs, include_novertgraph=True)
            ret3 = [
                GG for G in hairy for GG in self._hairy_to_w_edge_1eps1w(G)]
            # graphs with omega-tadpole, and valency of the eps vertex 1, not counting the tadpole
            hairy = self.get_hairy_graphs(
                self.n_vertices, self.n_loops-1, self.n_hairs, include_novertgraph=True)
            ret4 = [GG for G in hairy for GG in self._hairy_to_w_edge_1eps1w(
                G, addtadpole=True)]
            # graphs w/o omega tadpole, and valency of the eps vertex equals 0
            hairy = self.get_hairy_graphs(
                self.n_vertices, self.n_loops, self.n_hairs+1, include_novertgraph=True)
            ret5 = [GG for G in hairy for GG in self._hairy_to_w_hair_1w(G)]
            # graphs with omega-tadpole, and valency of the eps vertex 0, not counting the tadpole
            hairy = self.get_hairy_graphs(
                self.n_vertices, self.n_loops-1, self.n_hairs+1, include_novertgraph=True)
            ret6 = [GG for G in hairy for GG in self._hairy_to_w_hair_1w(
                G, addtadpole=True) for G in hairy]
            # missing 1w 2 eps graphs
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops, self.n_hairs+2)
            ret7 = [GG for G in hairy for GG in self._hairy_to_w_hairhair_1w2eps(
                G, addtadpole=False) for G in hairy]
            hairy = self.get_hairy_graphs(
                self.n_vertices+1, self.n_loops-1, self.n_hairs+2)
            ret8 = [GG for G in hairy for GG in self._hairy_to_w_hairhair_1w2eps(
                G, addtadpole=True) for G in hairy]
            # missing 1w 1eps graphs
            hairy = self.get_hairy_graphs(
                self.n_vertices, self.n_loops, self.n_hairs+2)
            ret9 = [GG for G in hairy for GG in self._hairy_to_w_hairhair_1w1eps(
                G, addtadpole=False) for G in hairy]
            hairy = self.get_hairy_graphs(
                self.n_vertices, self.n_loops-1, self.n_hairs+2)
            ret10 = [GG for G in hairy for GG in self._hairy_to_w_hairhair_1w1eps(
                G, addtadpole=True) for G in hairy]

        # Produce all permutations of the hairs
        all_perm = [list(range(self.n_vertices+2)) + list(p)
                    for p in itertools.permutations(range(self.n_vertices+2, self.n_vertices+self.n_hairs+2))]
        # return [G.relabel(p, inplace=False) for p in all_perm ]

        return [G.relabel(p, inplace=False) for p in all_perm for G in ret1 + ret2 + ret3 + ret4 + ret5 + ret6 + ret7 + ret8+ret9+ret10]

    def perm_sign(self, G, p):
        # The sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the hair-vertices.
        sgn = self.ogvs.perm_sign(G, p)

        # Compute the extra contribution from hairs.
        # if self.even_hairs == self.even_edges:
        #     hairs = p[self.n_vertices:]
        #     if len(hairs) != 0:
        #         sgn *= Shared.Perm.shifted(hairs).signature()
        return sgn

    def _bip_to_ordinary(self, G, nvertices, nedges, nhairs):
        # Translates bipartite into ordinary graph by replacing a bivalent vertex of colour 2 with an edge.
        for v in range(nvertices, nvertices + nhairs + nedges):
            neighbors = G.neighbors(v)
            n_l = len(neighbors)
            if n_l == 1:  # hair
                continue
            elif n_l == 2:  # edge
                G.add_edge(neighbors)
                G.delete_vertex(v)
            else:
                raise ValueError(
                    '%s: Vertices of second colour should have 1 or 2 neighbours' % str(self))
        G.relabel(range(G.order()))
        return G

    def _hairy_to_w_edge(self, G, addtadpole=False):
        # Translates hairy graph into w graph by taking vertices adjaent to edge as extra ones
        # Assumes we have n_vertices + 2 vertices in the graph
        for i in range(self.n_vertices + 2):
            for j in range(self.n_vertices + 2):
                if G.has_edge(i, j) and G.degree(i) == self.n_ws + 1:
                    # found valid combination: make our two vertices the last two int vertices
                    GG = G.relabel(dict_to_list(
                        {i: self.n_vertices + 1, self.n_vertices+1: i}, G.order()), inplace=False)
                    newj = i if j == self.n_vertices+1 else j
                    GG.relabel(dict_to_list(
                        {newj: self.n_vertices, self.n_vertices: newj}, GG.order()))
                    # j: self.n_vertices, self.n_vertices : j,
                    GG.delete_edge(self.n_vertices, self.n_vertices + 1)
                    yield GG

    def _hairy_to_w_edge_1eps(self, G, addtadpole=False):
        # Translates hairy graph into w graph by taking one vertex adjacent to edge as extra one...
        # with eps a new vertex of valency one
        # Assumes we have n_vertices + 1 vertices in the graph
        for i in range(self.n_vertices + 1):
            for j in range(self.n_vertices + self.n_hairs + 1):
                if G.has_edge(i, j) and G.degree(i) == self.n_ws + 1:
                    # found valid combination, relabel
                    GG = G.relabel(dict_to_list(
                        {i: self.n_vertices + 1, self.n_vertices+1: self.n_vertices+self.n_hairs+1}, G.order()), inplace=False)
                    GG.add_vertex(i)
                    GG.relabel(dict_to_list(
                        {i: self.n_vertices, self.n_vertices: i}, GG.order()))
                    newj = j
                    if j == self.n_vertices:
                        newj = i
                    elif j == self.n_vertices+1:
                        newj = self.n_vertices+self.n_hairs+1

                    GG.delete_edge(self.n_vertices+1, newj)
                    GG.add_edge(self.n_vertices, newj)
                    yield GG

    def _hairy_to_w_edge_1w(self, G, addtadpole=False):
        # Translates hairy graph into w graph by taking one vertex adjacent to edge as extra one...
        # Pendant to _hairy_to_w_edge_1eps
        for i in range(self.n_vertices + 1):
            for j in range(self.n_vertices + self.n_hairs + 1):
                if G.has_edge(i, j):
                    # found valid combination, relabel
                    GG = G.relabel(dict_to_list(
                        {i: self.n_vertices, self.n_vertices: i, self.n_vertices+1: self.n_vertices+self.n_hairs+1}, G.order()), inplace=False)
                    GG.add_vertex(self.n_vertices+1)
                    newj = j
                    if j == self.n_vertices:
                        newj = i
                    elif j == self.n_vertices+1:
                        newj = self.n_vertices+self.n_hairs+1

                    GG.delete_edge(self.n_vertices, newj)
                    GG.add_edge(self.n_vertices+1, newj)
                    yield GG

    def _hairy_to_w_edge_1eps1w(self, G, addtadpole=False):
        # Translates hairy graph into w graph by taking edge into two special hairs
        # with eps a new vertex of valency one
        # Assumes we have n_vertices vertices in the graph
        if self.n_ws != 1:
            return
        for i in range(self.n_vertices+self.n_hairs):
            for j in range(self.n_vertices+self.n_hairs):
                if G.has_edge(i, j):
                    # found valid combination, temporarily put eps and w at the end
                    GG = copy(G)
                    GG.add_vertex(self.n_vertices+self.n_hairs)
                    GG.add_vertex(self.n_vertices+self.n_hairs+1)
                    GG.delete_edge(i, j)
                    GG.add_edge(self.n_vertices+self.n_hairs, i)
                    GG.add_edge(self.n_vertices+self.n_hairs+1, j)
                    # relabel to move the vertices up front agai
                    GG.relabel({self.n_vertices: self.n_vertices+self.n_hairs, self.n_vertices+1: self.n_vertices+self.n_hairs+1,
                                self.n_vertices+self.n_hairs: self.n_vertices, self.n_vertices+self.n_hairs+1: self.n_vertices+1})
                    yield GG

    def _hairy_to_w_hair_2w(self, G, addtadpole=False):
        # Translates hairy graph into w graph
        # with eps a new vertex of valency zero
        # Assumes we have n_vertices + 1 vertices in the graph and n_hairs+1 hairs
        for i in range(self.n_vertices + 1):
            for j in range(self.n_vertices + 1, self.n_vertices+self.n_hairs+2):
                if G.has_edge(i, j) and G.degree(i) == self.n_ws+1:
                    # found valid combination, relabel such that i becomes w and j becomes eps
                    GG = copy(G)
                    GG.relabel(dict_to_list(
                        {i: self.n_vertices, j: self.n_vertices+1, self.n_vertices+1: j, self.n_vertices: i}, GG.order()))
                    GG.relabel(dict_to_list(
                        {self.n_vertices+1: self.n_vertices, self.n_vertices: self.n_vertices+1}, GG.order()))
                    GG.delete_edge(self.n_vertices, self.n_vertices + 1)
                    yield GG

    def _hairy_to_w_hairhair_2w(self, G, addtadpole=False):
        # Translates hairy graph into w graph
        # with eps a new vertex of valency one
        # Assumes we have n_vertices + 1 vertices in the graph and n_hairs+2 hairs
        # One hair points towards the new omega vertex, the other becomes eps
        # We can only consider graphs where both targets of the hairs are connected by an edge
        # because other graphs are covered by other routines
        for i in range(self.n_vertices + 1):
            for j in range(self.n_vertices + 1, self.n_vertices+self.n_hairs+3):
                if G.has_edge(i, j) and G.degree(i) == self.n_ws+1:
                    for k in range(self.n_vertices + 1, self.n_vertices+self.n_hairs+3):
                        if (k != j) and G.has_edge(i, G.neighbors(k)[0]):
                            # found valid combination, relabel such that i becomes w and k becomes eps
                            GG = copy(G)
                            GG.relabel(dict_to_list(
                                {i: self.n_vertices, k: self.n_vertices+1, self.n_vertices+1: k, self.n_vertices: i}, GG.order()))
                            GG.relabel(dict_to_list(
                                {self.n_vertices+1: self.n_vertices, self.n_vertices: self.n_vertices+1}, GG.order()))
                            #GG.delete_edge(self.n_vertices, self.n_vertices +1)
                            newj = k if j == self.n_vertices+1 else j
                            GG.delete_vertex(newj)
                            GG.relabel(
                                list(range(self.n_vertices+self.n_hairs+2)))
                            yield GG

    def _hairy_to_w_hairhair_1w2eps(self, G, addtadpole=False):
        # Translates hairy graph into w graph
        # with omega a vertex of valency one
        # Assumes we have n_vertices + 1 vertices in the graph and n_hairs+2 hairs
        # One hair points towards the new eps vertex, the other becomes w
        # This routine is the same as _hairy_to_w_hairhair_2w, just with eps and w interchanged
        for i in range(self.n_vertices + 1):
            for j in range(self.n_vertices + 1, self.n_vertices+self.n_hairs+3):
                if G.has_edge(i, j):
                    for k in range(self.n_vertices + 1, self.n_vertices+self.n_hairs+3):
                        if (k != j) and G.has_edge(i, G.neighbors(k)[0]):
                            # found valid combination, relabel such that i becomes w and k becomes eps
                            GG = copy(G)
                            GG.relabel(dict_to_list(
                                {i: self.n_vertices, k: self.n_vertices+1, self.n_vertices+1: k, self.n_vertices: i}, GG.order()))
                            GG.relabel(dict_to_list(
                                {self.n_vertices+1: self.n_vertices, self.n_vertices: self.n_vertices+1}, GG.order()))
                            #GG.delete_edge(self.n_vertices, self.n_vertices +1)
                            newj = k if j == self.n_vertices+1 else j
                            GG.delete_vertex(newj)
                            GG.relabel(
                                list(range(self.n_vertices+self.n_hairs+2)))
                            GG.relabel(
                                {self.n_vertices: self.n_vertices+1, self.n_vertices+1: self.n_vertices})
                            yield GG

    def _hairy_to_w_hairhair_1w1eps(self, G, addtadpole=False):
        # Translates hairy graph into w graph
        # with eps and omega of valency one
        # Assumes we have n_vertices vertices in the graph and n_hairs+2 hairs
        # One hair becomes the new omega vertex, the other becomes eps
        # We can only consider graphs where both targets of the hairs are connected by an edge
        # because other graphs are covered by other routines
        for i in range(self.n_vertices, self.n_vertices+self.n_hairs+2):
            for j in range(self.n_vertices, self.n_vertices+self.n_hairs+2):
                if G.has_edge(G.neighbors(i)[0], G.neighbors(j)[0]):
                    # found valid combination, relabel such that i becomes w and j becomes eps
                    GG = copy(G)
                    GG.relabel(dict_to_list(
                        {i: self.n_vertices+1, self.n_vertices+1: i}, GG.order()))
                    # GG.relabel(dict_to_list({self.n_vertices+1: self.n_vertices, self.n_vertices: self.n_vertices+1}, GG.order()))
                    #GG.delete_edge(self.n_vertices, self.n_vertices +1)
                    newj = i if j == self.n_vertices+1 else j
                    GG.relabel(dict_to_list(
                        {newj: self.n_vertices, self.n_vertices: newj}, GG.order()))
                    yield GG

    def _hairy_to_w_hair_1w(self, G, addtadpole=False):
        # Translates hairy graph into w graph
        # with eps a new vertex of valency zero
        # Assumes we have n_vertices vertices in the graph and n_hairs+1 hairs
        for j in range(self.n_vertices, self.n_vertices+self.n_hairs+1):
            # found valid combination, relabel such that i becomes w and j becomes eps
            GG = G.relabel({self.n_vertices: self.n_vertices +
                           self.n_hairs+1}, inplace=False)
            newj = self.n_vertices+self.n_hairs+1 if j == self.n_vertices else j
            GG.relabel({self.n_vertices+1: newj, newj: self.n_vertices+1})
            GG.add_vertex(self.n_vertices)
            yield GG


class WHairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of hairy graph vector spaces with specified number of omega hairs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_range (range): Range for the number of hairs.
        - w_range (range): number of omega hairs
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_range, w_range):
        """Initialize the sum vector space.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param h_range: Range for the number of hairs.
        :type h_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs: True for even hairs, False for odd hairs.
        :type even_hairs: bool
        """
        self.v_range = v_range
        self.l_range = l_range
        self.h_range = h_range
        self.w_range = w_range
        self.sub_type = ""

        vs_list = [WHairyGraphVS(v, l, h, w) for
                   (v, l, h, w) in itertools.product(self.v_range, self.l_range, self.h_range, self.w_range)]
        super().__init__(vs_list)

    def get_type(self):
        return 'whairy graphs'

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('hairs', self.h_range), ('ws', self.w_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s" % graph_type
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


# ------- Operators --------
class ContractEdgesGO(GraphOperator.GraphOperator):
    """Contract edges graph operator.

    Operate on a w-hairy graph by contracting an edge not connected to a hair vertex and unifying the two adjacent vertices.

    Attributes:
        - sub_type (str): Graphs sub type of the domain.
    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: HairyGraphVS
        :param target: Target vector space of the operator.
        :type target: HairyGraphVS
        """
        self.sub_type = domain.sub_type
        super().__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding contract edges graph operator.

        The contract edges operator reduces the number of vertices by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: HairyGraphVS
        :param target: Potential target vector space of the operator.
        :type target: HairyGraphVS
        :return: True if domain and target match to generate a corresponding contract edges graph operator.
        :rtype: bool
        """
        return domain.n_vertices == target.n_vertices + 1 and domain.n_loops == target.n_loops \
            and domain.n_hairs == target.n_hairs and domain.n_ws == target.n_ws

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, n_ws):
        """Return a contract edges graph operator.

        :param n_vertices: Number of vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param n_hairs: Number of hairs.
        :type n_hairs: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs: True for even hairs, False for odd hairs.
        :type even_hairs: bool
        :return: Contract edges graph operator based on the specified domain vector space.
        :rtype: ContractEdgesGO
        """
        domain = WHairyGraphVS(n_vertices, n_loops, n_hairs, n_ws)
        target = WHairyGraphVS(n_vertices - 1, n_loops, n_hairs, n_ws)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        if not self.is_valid():
            return 0
        try:
            (domain_dim, dimtarget_dim) = (
                self.domain.get_dimension(), self.target.get_dimension())
        except StoreLoad.FileNotFoundError:
            return 0
        if domain_dim == 0 or dimtarget_dim == 0:
            return 0
        return self.domain.n_edges * math.log(self.target.get_dimension(), 2)

    def get_type(self):
        return 'contract edges'

    def operate_on(self, G):
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image = []
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # only edges not connected to a numbered hair-vertex can be contracted
            if u >= self.domain.n_vertices+2 or v >= self.domain.n_vertices+2:
                continue

            # ensure u<v (this should be always true anyway actually)
            if u > v:
                u, v = v, u

            sgn = 1 if i % 2 == 0 else -1
            previous_size = G.size()
            previous_has_tadpole = (
                previous_size - self.domain.n_vertices - self.domain.n_hairs < self.domain.n_loops)
            sgn *= -1 if previous_has_tadpole else 1
            G1 = copy(G)
            # label all edges to determine sign later
            Shared.enumerate_edges(G1)

            # we always delete the lower index vertex. This ensures that the extra vertices are never deleted
            if v <= self.domain.n_vertices:
                G1.merge_vertices([v, u])
                if (previous_size - G1.size()) != 1:
                    continue
                G1.relabel(range(self.domain.n_vertices+1 +
                           self.domain.n_hairs), inplace=True)
                # find edge permutation sign
                sgn *= Shared.shifted_edge_perm_sign2(G1)
                image.append((G1, sgn))
            elif v == self.domain.n_vertices+1:
                # the second vertex is now omega, so we need to merge the vertex with the eps vertex
                # after reonnecting one of the edges to omega
                # we assume that u != eps, because this is forbidden in our graphs
                G1.delete_edge(u, v)
                # special care must be taken since a tadpole could be created at eps
                # and this is true iff there is an edge u-eps
                eps = self.domain.n_vertices
                new_has_tadpole = G1.has_edge(u, eps)
                # double tadpole => zero
                if new_has_tadpole and previous_has_tadpole:
                    continue
                if new_has_tadpole:
                    # remove the edge and compute the appropriate sign
                    k = G1.edge_label(u, eps)
                    G1.delete_edge(u, eps)
                    sgn *= 1 if ((k % 2 == 0) == (k < i)) else -1
                # loop over other neighbors w to be connected to omega
                for w in G1.neighbors(u):
                    G2 = copy(G1)
                    sgn2 = sgn
                    # reconnect the w-v-edge to omega (i.e., to v)
                    old_label = G2.edge_label(u, w)
                    G2.delete_edge(u, w)
                    G2.add_edge(w, v, old_label)

                    # now merge u and eps
                    G2.merge_vertices([eps, u])
                    # in case we have too few edges some double edges have been created => zero
                    if (previous_size - G2.size()) != (2 if new_has_tadpole else 1):
                        continue
                    G2.relabel(range(self.domain.n_vertices+1 +
                               self.domain.n_hairs), inplace=True)
                    # find edge permutation sign
                    sgn2 *= Shared.shifted_edge_perm_sign2(G2)
                    image.append((G2, sgn2))

        return image


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: HairyGraphSumVS
        """
        super().__init__(sum_vector_space,
                                             ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

# ------- Graph Complex --------


class WHairyGC(GraphComplex.GraphComplex):
    """Graph complex for hairy graphs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_range (range): Range for the number of hairs.
        - w_range (range): Range for the number of omega hairs
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_range, w_range, differentials):
        """Initialize the graph complex.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param h_range: Range for the number of hairs.
        :type  h_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs: True for even hairs, False for odd hairs.
        :type even_hairs: bool
        :param differentials: List of differentials. Options: 'contract', 'et1h'.
        :type differentials: list(str)
        """
        self.v_range = v_range
        self.l_range = l_range
        self.h_range = h_range
        self.w_range = w_range
        self.sub_type = ""

        sum_vector_space = WHairyGraphSumVS(
            self.v_range, self.l_range, self.h_range, self.w_range)
        differential_list = []
        if not set(differentials).issubset(['contract']):
            raise ValueError(
                "Differentials for hairy graph complex: 'contract'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)
        super().__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))

    def print_dim_and_eulerchar(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    ds = [WHairyGraphVS(v, l, h, w).get_dimension()
                          for v in self.v_range]
                    eul = sum([(1 if j % 2 == 0 else -1) *
                              d for j, d in enumerate(ds)])
                    print("Dimensions (w,h,l) ", w,
                          h, l, ":", ds, "Euler", eul)

    def print_cohomology_dim(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    cohomdict = {}
                    for v in self.v_range:
                        D1 = ContractEdgesGO.generate_operator(v, l, h, w)
                        D2 = ContractEdgesGO.generate_operator(v+1, l, h, w)
                        try:
                            d = WHairyGraphVS(v, l, h, w).get_dimension()
                            r1 = D1.get_matrix_rank()
                            r2 = D2.get_matrix_rank()
                            cohomdict[v] = d-r1-r2
                        except:
                            pass

                    print("Cohomology Dimensions (w,h,l) ",
                          w, h, l, ":", cohomdict)
