### This file contains experiments on a tentative cyclic version of the graph complex


"""Graph complexes based on ordinary (simple, without multiple edges) graphs.
Implemented Differentials: Contract edges, delete edges.
The generators only produce 1-vertex irreducible graphs."""


__all__ = ['graph_type', 'sub_types', 'OrdinaryGVS', 'OrdinaryGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
           'DeleteEdgesGO', 'DeleteEdgesD', 'OrdinaryGC']

import itertools
from sage.all import *
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import NautyInterface
import Parameters
import GCDimensions
import HairyGraphComplex
import OrdinaryGraphComplex


graph_type = "gograph"



# ------- Graph Vector Space --------
class GOneVS(GraphVectorSpace.GraphVectorSpace):
    """Ordinary graph vector space, with one marked vertex, i.e., Graphs_2(1).
    The marked vertex is the first in the ordering.


    Attributes:
        - n_vertices (int): Number of internal vertices.
        - n_loops (int): Number of loops.
        - n_edges (int): Number of edges.

    """

    def __init__(self, n_vertices, n_loops):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_edges = self.n_loops + self.n_vertices
        super(GOneVS, self).__init__()

    def get_type(self):
        return '%s graphs' % (graph_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops

    def __hash__(self):
        return hash("gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return [[0], list(range(1,self.n_vertices+1))]

    def is_valid(self):
        # internal Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return (3 * self.n_vertices <= 2 * self.n_edges) and self.n_vertices >= 0 and self.n_loops >= 0 \
            and self.n_edges <= self.n_vertices * (self.n_vertices + 1) / 2

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return GCDimensions.get_ordinary_dim_estimate(self.n_vertices+1, self.n_loops)
        #return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        if self.n_vertices == 0 and self.n_loops == 0:
            # Output the empty graph
            G = Graph(1)
            return [ G ]

        # for ext_valence in range(1,self.n_vertices+1):
        #     # take hairy graph complex
        #     # ideally, hairy gc should have been created including non-1vi vertices
        #     HGC = HairyGraphComplex.HairyGraphVS(self.n_vertices, self.n_loops - ext_valence+1, ext_valence, False, True)
        #     for G in HGC.get_basis():
        #         yield G

        for G in NautyInterface.list_simple_graphs_1(self.n_vertices+1, self.n_edges, onlyonevi=False):
            for j in range(self.n_vertices+1):
                p = [ (i+j) % (self.n_vertices+1) for i in range(self.n_vertices+1) ]
                GG = G.relabel(p, inplace=False)
                # check if valence conditions satisfied
                if all( GG.degree(v+1) >=3 for v in range(self.n_vertices) ):
                    yield GG

    def perm_sign(self, G, p):
        # The sign is (induced sign of the edge permutation)
        # We assume the edges are always lexicographically ordered
        # For the computation we use that G.edges() returns the edges in lex ordering
        # We first label the edges on a copy of G lexicographically
        G1 = copy(G)
        Shared.enumerate_edges(G1)
        # We permute the graph, and read of the new labels
        G1.relabel(p, inplace=True)
        return Shared.Perm([j for (u, v, j) in G1.edges()]).signature()




# ------- Operators --------
class ReconnectEdgesGO(GraphOperator.GraphOperator):
    """Reconnect edges graph operator.

    Operates on an ordinary graph with one marked vertex by 
    reconnecting an arbitrary subset of the half-edges to the special vertex.
    Afterwards makes the special vertex normal, so that the image is an ordinary graph.

    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: OrdinaryGVS
        :param target: Target vector space of the operator.
        :type target: OrdinaryGVS
        """
        if not ReconnectEdgesGO.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        super(ReconnectEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding contract edges graph operator.

        The contract edges operator reduces the number of vertices by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: OrdinaryGVS
        :param target: Potential target vector space of the operator.
        :type target: OrdinaryGVS
        :return: True if domain and target match to generate a corresponding contract edges graph operator.
        :rtype: bool
        """
        return domain.n_vertices + 1 == target.n_vertices and domain.n_loops == target.n_loops \
            and target.even_edges == False

    @classmethod
    def generate_operator(cls, n_vertices, n_loops):
        """Return a contract edge graph operator.

        :param n_vertices: Number of (internal) vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :return: Contract edges graph operator based on the specified domain vector space.
        :rtype:ContractEdgesGO
        """
        domain = GOneVS(n_vertices, n_loops)
        target = OrdinaryGraphComplex.OrdinaryGVS(n_vertices + 1, n_loops, False)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "reconnectF%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "reconnectF%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_matrix_file_path(self):
        s = "reconnectF%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ref_rank_file_path(self):
        s = "reconnectF%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        if not self.is_valid():
            return 0
        return self.domain.n_edges # think: 2**edges

    def get_type(self):
        return 'reconnect edges'


    def operate_on(self, G):
        # Operates on the graph G by reconnecting half-edges to the special vertex 0
        image = []

        def reconnect_rec(GG, e_list, sgn, image):
            # recursive subroutine that reconnects the edges in e_list in all possible ways
            if len(e_list) == 0:
                # output graph
                sgn2 = sgn * Shared.shifted_edge_perm_sign(GG)
                image.append((G1, sgn2))
            else:
                (u,v) = e_list[0]
                ### three ways: connect none, first, or second half-edge
                # no reconnect
                reconnect_rec(GG, e_list[1:], sgn, image)
                sgn *= -1
                # reconnect first half-edge 
                if u != 0 and v != 0 and GG.degree(u) >= 4:
                    G1 = copy(GG)
                    lbl = G1.edge_label(u,v)
                    G1.delete_edge(u,v)
                    G1.add_edge(0,v,label=lbl)
                    reconnect_rec(G1, e_list[1:], sgn, image)
                # reconnect second half-edge 
                if u != 0 and v != 0 and GG.degree(v) >= 4:
                    G2 = copy(GG)
                    lbl = G2.edge_label(u,v)
                    G2.delete_edge(u,v)
                    G2.add_edge(0,u,label=lbl)
                    reconnect_rec(G2, e_list[1:], sgn, image)

        reconnect_rec(G, G.edges(labels=False), 1, image)
        # for (i, e) in enumerate(G.edges(labels=False)):
        #     (u, v) = e
        #     # print("contract", u, v)
        #     pp = Shared.permute_to_left(
        #         (u, v), range(0, self.domain.n_vertices))
        #     sgn = self.domain.perm_sign(G, pp)
        #     G1 = copy(G)
        #     G1.relabel(pp, inplace=True)
        #     Shared.enumerate_edges(G1)
        #     previous_size = G1.size()
        #     G1.merge_vertices([0, 1])
        #     if (previous_size - G1.size()) != 1:
        #         continue
        #     # print(sgn)
        #     G1.relabel(list(range(0, G1.order())), inplace=True)
        #     if not self.domain.even_edges:
        #         # p = [j for (a, b, j) in G1.edges()]
        #         # sgn *= Permutation(p).signature()
        #         sgn *= Shared.shifted_edge_perm_sign(G1)
        #     else:
        #         sgn *= -1  # TODO overall sign for even edges
        #     image.append((G1, sgn))
        return image


