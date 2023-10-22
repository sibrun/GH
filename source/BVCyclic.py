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


# whether to allow disconnect graphs in the definition of GOneVS
allow_disconnected = True
alternate_generation = False

graph_type = "gograph" + ("_d" if allow_disconnected else "")+ ("_a" if alternate_generation else "")




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
        if alternate_generation:
            return self.get_generating_graphs2()
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

        if allow_disconnected:
            # add all graphs consisting of a single external vertex and a disconnected internal component
            for G in NautyInterface.list_simple_graphs(self.n_vertices, self.n_edges, onlyonevi=False):
                G.add_vertex()
                # make new vertex the 0-th
                p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
                G.relabel(p, inplace=True)
                yield G

        for G in NautyInterface.list_simple_graphs_1(self.n_vertices+1, self.n_edges, onlyonevi=False):
            for j in range(self.n_vertices+1):
                p = [ (i+j) % (self.n_vertices+1) for i in range(self.n_vertices+1) ]
                GG = G.relabel(p, inplace=False)
                # check if valence conditions satisfied
                if all( GG.degree(v+1) >=3 for v in range(self.n_vertices) ):
                    yield GG

    def get_generating_graphs2(self):
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

        if allow_disconnected:
            # add all graphs consisting of a single external vertex and a disconnected internal component
            for G in NautyInterface.list_simple_graphs(self.n_vertices, self.n_edges, onlyonevi=False):
                G.add_vertex()
                # make new vertex the 0-th
                p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
                G.relabel(p, inplace=True)
                yield G

        # add all graphs with >=trivalent external vertex
        for G in NautyInterface.list_simple_graphs(self.n_vertices+1, self.n_edges, onlyonevi=False):
            for j in range(self.n_vertices+1):
                p = [ (i+j) % (self.n_vertices+1) for i in range(self.n_vertices+1) ]
                GG = G.relabel(p, inplace=False)
                yield GG
        
        # add graphs with uni- or bivalent external vertex using hairy graphs
        hvs1 = HairyGraphComplex.HairyGraphVS(self.n_vertices, self.n_loops, 1, False, True)
        for G in hvs1.get_basis():
            # relabel so that hair is first vertex
            p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
            G.relabel(p, inplace=True)
            yield G
        hvs2 =  HairyGraphComplex.HairyGraphVS(self.n_vertices, self.n_loops-1, 2, False, True)
        for G in hvs2.get_basis():
            # merge hairs
            G.merge_vertices([self.n_vertices, self.n_vertices+1])
            G.relabel(list(range(0, self.n_vertices+1)), inplace=True)
            # relabel so that hairs are first vertex
            p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
            G.relabel(p, inplace=True)
            yield G


        # add all graphs with univalent external vertex, connected to a >=4-valent vertex
        # ...and all graphs with bivalent external vertex
        # for G in NautyInterface.list_simple_graphs(self.n_vertices, self.n_edges-1, onlyonevi=False):
        #     G.add_vertex()
        #     # make new vertex the 0-th
        #     p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
        #     G.relabel(p, inplace=True)
        #     # add an edge
        #     for j in range(self.n_vertices):
        #         GG = copy(G)
        #         GG.add_edge(0,j+1)
        #         yield GG
        #     # or put 0 in the middle of an edge
        #     for u,v in G.edges(labels=False):
        #         GG = copy(G)
        #         GG.delete_edge( (u,v) )
        #         GG.add_edge(0,u)
        #         GG.add_edge(0,v)
        #         yield GG

        # add all graphs with univalent external vertex, connected to a 3-valent vertex
        # for G in NautyInterface.list_simple_graphs(self.n_vertices-1, self.n_edges-2, onlyonevi=False):
        #     G.add_vertex()  # new 0 vertex
        #     G.add_vertex()  # new neighbor of 0 -> vertex 1
        #     # make new vertex the 0-th
        #     p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
        #     G.relabel(p, inplace=True)
        #     # iterate over all edges and connect 0 to the middle of the edge
        #     for u,v in G.edges(labels=False):
        #         GG = copy(G)
        #         GG.delete_edge( (u,v) )
        #         GG.add_edge(0,1)
        #         GG.add_edge(u,1)
        #         GG.add_edge(v,1)
        #         yield GG
        # # ... plus cases where a bivalent vertex is added parallel to an existing edge
        # for G in NautyInterface.list_simple_graphs(self.n_vertices, self.n_edges-2, onlyonevi=False):
        #     G.add_vertex()  # new 0 vertex
        #     # make new vertex the 0-th
        #     p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
        #     G.relabel(p, inplace=True)
        #     # iterate over all edges and connect 0 to the middle of the edge
        #     for u,v in G.edges(labels=False):
        #         GG = copy(G)
        #         GG.add_edge(u,0)
        #         GG.add_edge(v,0)
        #         yield GG

        # # ...plus cases where we put the new double-edge parallel to an existing one
        # for G in NautyInterface.list_simple_graphs(self.n_vertices-1, self.n_edges-3, onlyonevi=False):
        #     G.add_vertex()  # new 0 vertex
        #     G.add_vertex()  # new neighbor of 0 -> vertex 1
        #     # make new vertex the 0-th
        #     p = [self.n_vertices - j for j in range(self.n_vertices+1) ]
        #     G.relabel(p, inplace=True)
        #     # iterate over all edges and connect 0 to the middle of a parallel edge
        #     for u,v in G.edges(labels=False):
        #         GG = copy(G)
        #         GG.add_edge(0,1)
        #         GG.add_edge(u,1)
        #         GG.add_edge(v,1)
        #         yield GG


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


class GOneVS3V(GraphVectorSpace.GraphVectorSpace):
    """Ordinary graph vector space, with one marked vertex, i.e., Graphs_2(1).
    The marked vertex is the first in the ordering.
    In comparison to GOneVS we require here that the exernal vertex has valence at least 3


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
        self.gonevs = GOneVS(n_vertices, n_loops)
        super(GOneVS3V, self).__init__()

    def get_type(self):
        return '%s graphs' % (graph_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops

    def __hash__(self):
        return hash("gra3v%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra3v%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "gra3v%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return [[0], list(range(1,self.n_vertices+1))]

    def is_valid(self):
        # internal Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return (3 * self.n_vertices +3 <= 2 * self.n_edges) and self.n_vertices >= 0 and self.n_loops >= 0 \
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
    
        for G in NautyInterface.list_simple_graphs(self.n_vertices+1, self.n_edges, onlyonevi=False):
            for j in range(self.n_vertices+1):
                p = [ (i+j) % (self.n_vertices+1) for i in range(self.n_vertices+1) ]
                GG = G.relabel(p, inplace=False)
                yield GG

    def perm_sign(self, G, p):
        # The sign is (induced sign of the edge permutation)
        # We assume the edges are always lexicographically ordered
        # For the computation we use that G.edges() returns the edges in lex ordering
        # We first label the edges on a copy of G lexicographically
        return self.gonevs.perm_sign(G,p)





# ------- Operators --------
class ReconnectEdgesGO(GraphOperator.GraphOperator):
    """Reconnect edges graph operator.

    Operates on an ordinary graph with one marked vertex by 
    reconnecting an arbitrary subset of the half-edges to the special vertex.
    Afterwards makes the special vertex normal, so that the image is an ordinary graph.

    The image is multiplied by (-1)**(degree of vertex 0).
    Also, the identity matrix is implicitly added. 

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
        # start with identity operator
        image.append( (G, 1) )

        def reconnect_rec(GG, e_list, sgn, image):
            # recursive subroutine that reconnects the edges in e_list in all possible ways
            if len(e_list) == 0:
                # output graph
                sgn2 = sgn * Shared.shifted_edge_perm_sign(GG)
                image.append((GG, sgn2))
            else:
                (u,v) = e_list[0]
                ### three ways: connect none, first, or second half-edge
                # no reconnect
                reconnect_rec(GG, e_list[1:], sgn, image)
                sgn *= -1
                # reconnect first half-edge 
                if u != 0 and v != 0 and GG.degree(u) >= 4 and not GG.has_edge(0,v):
                    G1 = copy(GG)
                    lbl = G1.edge_label(u,v)
                    G1.delete_edge(u,v)
                    G1.add_edge(0,v,label=lbl)
                    reconnect_rec(G1, e_list[1:], sgn, image)
                # reconnect second half-edge 
                if u != 0 and v != 0 and GG.degree(v) >= 4 and not GG.has_edge(0,u):
                    G2 = copy(GG)
                    lbl = G2.edge_label(u,v)
                    G2.delete_edge(u,v)
                    G2.add_edge(0,u,label=lbl)
                    reconnect_rec(G2, e_list[1:], sgn, image)

        Gp = copy(G)
        Shared.enumerate_edges(Gp)
        reconnect_rec(Gp, Gp.edges(labels=False), (-1) ** (Gp.degree(0)), image)
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

class ReconnectEdgesGO3V(GraphOperator.GraphOperator):
    """Reconnect edges graph operator.

    Operates on an ordinary graph with one marked vertex by 
    reconnecting an arbitrary subset of the half-edges to the special vertex.
    Afterwards makes the special vertex normal, so that the image is an ordinary graph.

    The image is multiplied by (-1)**(degree of vertex 0).
    Also, the identity matrix is implicitly added. 

    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: OrdinaryGVS
        :param target: Target vector space of the operator.
        :type target: OrdinaryGVS
        """
        if not ReconnectEdgesGO3V.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.reconnect_op = ReconnectEdgesGO.generate_operator(domain.n_vertices, domain.n_loops)
        super(ReconnectEdgesGO3V, self).__init__(domain, target)

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
        domain = GOneVS3V(n_vertices, n_loops)
        target = OrdinaryGraphComplex.OrdinaryGVS(n_vertices + 1, n_loops, False)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "reconnect3vF%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "reconnect3vF%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_matrix_file_path(self):
        s = "reconnect3vF%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ref_rank_file_path(self):
        s = "reconnect3vF%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        if not self.is_valid():
            return 0
        return self.domain.n_edges # think: 2**edges

    def get_type(self):
        return 'reconnect edges 3v'


    def operate_on(self, G):
        return self.reconnect_op.operate_on(G)


class AddVReconnectEdgesGO(GraphOperator.GraphOperator):
    """Reconnect edges graph operator of Marko Zivkovic.

    Operates on an ordinary graph by adding a vertex and then 
    reconnecting an arbitrary subset of the half-edges to the new vertex.
   
    The image is multiplied by (-1)**(degree of vertex 0).

    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: OrdinaryGVS
        :param target: Target vector space of the operator.
        :type target: OrdinaryGVS
        """
        if not AddVReconnectEdgesGO.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for addvreconnect  edges operator")
        
        self.reconnect_op = ReconnectEdgesGO.generate_operator(domain.n_vertices, domain.n_loops)
        super(AddVReconnectEdgesGO, self).__init__(domain, target)

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
        return domain.n_vertices + 1 == target.n_vertices and domain.n_loops == target.n_loops+1 \
            and target.even_edges == False and domain.even_edges == False

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
        domain = OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, False)
        target = OrdinaryGraphComplex.OrdinaryGVS(n_vertices + 1, n_loops-1, False)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "addvreconnectF%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "addvreconnectF%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_matrix_file_path(self):
        s = "addvreconnectF%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ref_rank_file_path(self):
        s = "addvreconnectF%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        if not self.is_valid():
            return 0
        return self.domain.n_edges # think: 2**edges

    def get_type(self):
        return 'add v and reconnect edges'


    def operate_on(self, G):
        # act by re-using ReconnectEdgesGO, on a graph built by adding one more external vertex
        # It is not necessary to subtract id, since that graph is not connected and automatically removed

        # add a new vertex that becomes the first
        G1 = copy(G)
        G1.add_vertex()
        p = [self.domain.n_vertices - j for j in range(self.domain.n_vertices+1) ]
        G1.relabel(p, inplace=True)
        return self.reconnect_op.operate_on(G1)




class GOneDegSlice(GraphVectorSpace.DegSlice):
    """Represents the sum of GOne and ordinary graph vector spaces used for the cohomology computation.
    n_vertices refers to the number of vertices in the ordinary GC part
    top_n=0 -> just ordinary GC
    top_n=1 -> also GOne present
    """

    def is_complete(self):
        for vs in self.vs_list:
            if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
                return False
        return True

    def __init__(self,  n_vertices, n_loops, top_n):
        """Initialize the degree slice.
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.top_n = top_n
        if top_n == 0:
            vs_list = [ OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, False)]
        elif top_n == 1:
            vs_list = [GOneVS(n_vertices-2, n_loops) , OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, False)]
        else:
            raise ValueError("GOneDegSlice init: top_n must be 0 or 1")
        
        super(GOneDegSlice, self).__init__(vs_list, n_vertices)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('topn', self.top_n)])

    def __eq__(self, other):
        return self.n_loops == other.n_loops \
            and self.n_vertices == other.n_vertices \
            and self.top_n == other.top_n

    def __str__(self):
        return ("GOneDegSlice_%s_%s_%s" % self.get_ordered_param_dict().get_value_tuple())

    def __hash__(self):
        return hash(str(self))

    def get_info_plot_path(self):
        s = "info_vertex_loop_top_degree_slice_deg_%d_%d_%d_%d_%s_%s" % (
            self.n_vertices, self.n_loops, self.top_n, graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


class ContractReconnectBiOM(GraphOperator.BiOperatorMatrix):
    def __init__(self, domain, target):
        """Bi operator matrix based on contract edges and reconnect edges.

        """
        self.domain = domain 
        self.target = target
        super(ContractReconnectBiOM, self).__init__(domain, target, OrdinaryGraphComplex.ContractEdgesGO,
                                                ReconnectEdgesGO)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target degree slices match to generate a corresponding bi operator matrix.

        The bi operator reduces the degree by one and increases the minimal number of hairs by one for both hair colors.

        :param domain: Potential domain vector space of the operator.
        :type domain: VertexLoopDegSlice
        :param target: Potential target vector space of the operator.
        :type target: VertexLoopDegSlice
        :return: True if domain and target match to generate a corresponding bi operator matrix.
        :rtype: bool
        """
        return domain.n_vertices == target.n_vertices + 1 \
            and domain.n_loops == target.n_loops

    def get_matrix_file_path(self):
        s = "bi_D_contract_reconnect_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "bi_D_contract_reconnect_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)
    
    @classmethod
    def generate_operator(cls, n_vertices, n_loops):
        domain = GOneDegSlice(n_vertices, n_loops,1)
        target = GOneDegSlice(n_vertices-1, n_loops,0)
        return cls(domain, target)


class GOneSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of forested graph vector spaces with specified edge parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range):
        self.l_range = l_range
        self.v_range = v_range
        # todo : misses top vs.... but seems ok in practice
        vs_list = [OrdinaryGraphComplex.OrdinaryGVS(v, l, False) for v in v_range for l in l_range] \
            + \
            [GOneVS(v-2, l) for v in v_range for l in l_range]

        super(GOneSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs' % (graph_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_top_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, s)


class ContractReconnectTopD(GraphOperator.Differential):
    """ Represents a collection of ContractUnmarkTopBiOM.
    This class is also used to compute cohomology.
    """

    def __init__(self, v_range, l_range):
        self.l_range = l_range
        self.v_range = v_range
        op_list = [ContractReconnectBiOM.generate_operator(v, l,)
                   for v in v_range
                   for l in l_range]
        sum_vs = GOneSumVS(v_range, l_range)
        super(ContractReconnectTopD, self).__init__(sum_vs, op_list)

    def get_type(self):
        return 'contract edges and reconnect edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_reconnect_edges_D_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_reconnect_edges_D_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        s = self.sum_vector_space
        return Shared.OrderedDict([('vertices', s.v_range), ('loops', s.l_range)])

    # def _get_single_cohomology(self, n_loops, n_marked, n_hairs):
    #     """
    #     Computes a single cohomology dimension.
    #     """
    #     opD = ContractUnmarkTopBiOM.generate_operator(
    #         n_loops, n_marked, n_hairs, self.even_edges)
    #     opDD = ContractUnmarkTopBiOM.generate_operator(
    #         n_loops, n_marked+1, n_hairs, self.even_edges)
    #     n_vertices = opD.domain.n_vertices

    #     opC = ContractEdgesGO.generate_operator(
    #         n_vertices, n_loops, n_marked, n_hairs, self.even_edges)

    #     try:
    #         # dimension of the kernel of the contraction
    #         dimV = opD.get_domain().get_dimension() - opC.get_matrix_rank()
    #     except StoreLoad.FileNotFoundError:
    #         logger.info(
    #             "Cannot compute cohomology: First build basis for %s " % str(opD.get_domain()))
    #         return None
    #     if dimV == 0:
    #         return '*'
    #     if opD.is_valid():
    #         try:
    #             rankD = opD.get_true_rank()
    #         except StoreLoad.FileNotFoundError:
    #             logger.info(
    #                 "Cannot compute cohomology: Matrix rank not calculated for %s " % str(opD))
    #             return None
    #     else:
    #         rankD = 0
    #     if opDD.is_valid():
    #         try:
    #             rankDD = opDD.get_true_rank()
    #         except StoreLoad.FileNotFoundError:
    #             logger.info(
    #                 "Cannot compute cohomology: Matrix rank not calculated for %s " % str(opDD))
    #             return None
    #     else:
    #         rankDD = 0
    #     cohomology_dim = dimV - rankD - rankDD
    #     if cohomology_dim < 0:
    #         raise ValueError("Negative cohomology dimension for %s (%d - %d - %d)" %
    #                          (str(opD.domain), dimV, rankD, rankDD))
    #         # logger.error("Negative cohomology dimension for %s" % str(opD.domain))
    #     return cohomology_dim

    # def get_cohomology_dim_dict(self):
    #     ms = list(self.m_range)[:-1]
    #     return {(l, m, h): self._get_single_cohomology(l, m, h)
    #             for l in self.l_range
    #             for h in self.h_range
    #             for m in ms}

        # return super().get_cohomology_dim_dict()

    # def build_basis(self, **kwargs):
    #     # self.sum_vector_space.compute_all_pregraphs(**kwargs)
    #     self.sum_vector_space.build_basis(**kwargs)

    # def compute_rank(self, sage=None, linbox=None, rheinfall=None, sort_key='size', ignore_existing_files=False, n_jobs=1, info_tracker=False):
    #     # compute ranks of contractto operators
    #     super().compute_rank(sage, linbox, rheinfall, sort_key,
    #                          ignore_existing_files, n_jobs, info_tracker)
    #     # compute ranks of contract operators that are also necessary to have
    #     print("Computing contract operator ranks...")
    #     coplist = [ContractEdgesGO.generate_operator(2*l-2+h, l, m, h, self.even_edges)
    #                for l in self.l_range
    #                for m in self.m_range
    #                for h in self.h_range]
    #     oc = GraphOperator.OperatorMatrixCollection(
    #         self.sum_vector_space, coplist)
    #     oc.compute_rank(sage, linbox, rheinfall, sort_key,
    #                     ignore_existing_files, n_jobs, info_tracker)

