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
from itertools import permutations
import OrdinaryGraphComplex

graph_type = "kneissler"

sub_types = {True: "even_edges", False: "odd_edges"}

# --- helpers for graph generation ---
def barrel_graph(k, p):
    # generates the barrel graph of 2k vertices, with p a permutation of the numbers 0,..,k-2
    G = Graph(2*k)
    # generate rims of barrel
    for j in range(k):
        G.add_edge(j, (j+1) % k)
        G.add_edge(k+j, k+(j+1)%k)
    # generate spokes
    G.add_edge(k-1, 2*k-1)
    for i,j in enumerate(p):
        G.add_edge(i, k+j)

    return G

def all_barrel_graphs(k):
    # generates all barrel graphs of 2k vertices
    for p in permutations(range(k-1)):
        yield barrel_graph(k, p)

def tbarrel_graph(k,p):
    # barrel graph with one 4-valent vertex (total 2k-1 vertices)
    G = barrel_graph(k, p)
    G.merge_vertices([k-2,k-1])
    G.relabel(list(range(0, G.order())), inplace=True)
    return G

def all_tbarrel_graphs(k):
    # generates all barrel graphs of 2k vertices
    for p in permutations(range(k-1)):
        yield tbarrel_graph(k, p)

def xtbarrel_graph(k,p):
    # handle barrel graph with one 4-valent vertex (total 2k-1 vertices)
    # p must be permutation of [0,..,k-2]
    G = Graph(2*k-1)
    # generate rims of barrel -- first rim is length k-1, second length k-1 
    for j in range(k-1):
        G.add_edge(j, (j+1) % (k-1))
    for j in range(0, k-1):
        G.add_edge(k+j, k+(j+1)%(k-1))
    # vertex k-1 is the "central" vertex
    # generate spokes
    G.add_edge(k-1, 2*k-2)
    G.add_edge(k-1, k-2)

    for i,j in enumerate(p):
        if j < k-2:
            G.add_edge(i, k+j)
        else: # j = k-2
            G.add_edge(i, k-1)

    return G

def all_xtbarrel_graphs(k):
    # generates all barrel graphs of 2k vertices
    for p in permutations(range(k-1)):
        if p[k-2] != k-2:
            yield xtbarrel_graph(k, p)


def triangle_graph(k,p):
    # trivalent graph, barrel with one central vertex, 2*k vertices
    # p must be permutation of [0,..,k-2]
    G = Graph(2*k)
    # generate rims of barrel -- first rim is length k, second length k-1 
    for j in range(k):
        G.add_edge(j, (j+1) % k)
    for j in range(k-1):
        G.add_edge(k+1+j, k+1+(j+1)%(k-1))
    # vertex k is the "central" vertex
    # generate spokes
    G.add_edge(k-1, k)
    G.add_edge(k, 2*k-1)

    for i,j in enumerate(p):
        if j < k-2:
            G.add_edge(i, k+1+j)
        else: # j = k-2
            G.add_edge(i, k)

    return G

def all_triangle_graphs(k):
    # generates all barrel graphs of 2k vertices
    for p in permutations(range(k-1)):
        yield triangle_graph(k, p)

def hgraph(k,p):
    # trivalent graph, barrel with two central vertices, 2*k vertices
    # p must be permutation of [0,..,k-2]
    G = Graph(2*k)
    # generate rims of barrel -- first rim is length k-1, second length k-1 
    for j in range(k-1):
        G.add_edge(j, (j+1) % (k-1))
    for j in range(k-1):
        G.add_edge(k+1+j, k+1+(j+1)%(k-1))
    # vertice k-1, k are the "central" vertex
    # generate spokes
    G.add_edge(k-2, k-1)
    G.add_edge(k-1, 2*k-1)
    G.add_edge(k-1, k)

    for i,j in enumerate(p):
        if i < k-2:
            G.add_edge(i, k+j)
        else: # i = k-2
            if j>0: # j=0 yields an invalid graph, that must be discarded
                G.add_edge(k, k+j)
    return G
def all_hgraph_graphs(k):
    # generates all barrel graphs of 2k vertices
    for p in permutations(range(k-1)):
        if p[k-2] > 0:
            yield hgraph(k, p)


# ------- Graph Vector Space --------
class KneisslerGVS(GraphVectorSpace.GraphVectorSpace):
    """Ordinary graph vector space.

    Sub vector space with specified number of vertices, loops and even or odd edges and at least trivalent vertices.
    No multiple edges. Only 1-vi graphs are produced by generators.

    Attributes:
        - n_vertices (int): Number of vertices.
        - n_loops (int): Number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - n_edges (int): Number of edges.
        - sub_type (str): Sub type of graphs.

    """

    def __init__(self, n_loops, type, even_edges):
        """Initialize the kneissler graph vector space.
        Possible types: 
        0 -- trivalent generators (barrel graphs)
        1 -- relation generators with one vertex of degree 4
        2 -- all graphs that can be produced by splitting a type 1 graph (trivalent H-graphs, triangle graphs or barrel graphs)
        3 -- the complement of type 0 within type 2.

        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param type: int: Type of graph.
        :type type: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        if type not in [0, 1, 2, 3]:
            raise ValueError("Type must be 0, 1, 2 or 3")
        self.k = n_loops - 1
        if type in [0,2,3]:
            self.n_vertices = 2 * self.k
        else:
            self.n_vertices = 2 * self.k - 1

        self.type = type
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(self.n_vertices, self.n_loops, self.even_edges)
        super().__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_loops == other.n_loops and self.type == other.type and self.even_edges == other.even_edges

    def __hash__(self):
        return hash("gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('loops', self.n_loops), ('type', self.type)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return self.type in [0, 1, 2, 3] and self.n_loops >= 0 

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return GCDimensions.get_ordinary_dim_estimate(self.n_vertices, self.n_loops)
        #return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        
        if self.type == 0:
            # trivalent generators
            for G in all_barrel_graphs(self.k):
                yield G
        elif self.type == 1:
            # relation generators
            for G in all_tbarrel_graphs(self.k):
                yield G
            for G in all_xtbarrel_graphs(self.k):
                yield G
        elif self.type == 2:
            # all graphs that can be produced by splitting a type 1 graph
            for G in all_barrel_graphs(self.k):
                yield G
            for G in all_triangle_graphs(self.k):
                yield G
            for G in all_hgraph_graphs(self.k):
                yield G
        elif self.type == 3:
            # the complement of type 0 within type 2
            # assume that types 0 and 2 have been computed
            V1 = KneisslerGVS(self.n_loops, 0, self.even_edges)
            V2 = KneisslerGVS(self.n_loops, 2, self.even_edges)
            b1 = V1.get_basis_g6()
            b2 = V2.get_basis_g6()
            # b1 is a subset of b2
            for g6 in set(b2) - set(b1):
                yield Graph(g6)


    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G, p)





# ------- Operators --------
class ContractEdgesGO(GraphOperator.GraphOperator):
    """Contract edges graph operator.

    Operates on an ordinary graph by contracting an edge and unifying the two adjacent vertices.

    Attributes:
        - sub_type (str): Graphs sub type of the domain.
    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: OrdinaryGVS
        :param target: Target vector space of the operator.
        :type target: OrdinaryGVS
        """
        if not ContractEdgesGO.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        self.ocontract = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(domain.n_vertices, domain.n_loops, domain.even_edges)
        super().__init__(domain, target)

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
        return domain.n_vertices - 1 == target.n_vertices and domain.n_loops == target.n_loops \
            and domain.even_edges == target.even_edges

    @classmethod
    def generate_operator(cls, n_loops, type, even_edges):
        """Return a contract edge graph operator.

        :param n_vertices: Number of vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :return: Contract edges graph operator based on the specified domain vector space.
        :rtype:ContractEdgesGO
        """
        if type not in [0, 2, 3]:
            raise ValueError("Type must be 0, 2 or 3")
        domain = KneisslerGVS(n_loops, type, even_edges)
        target = KneisslerGVS(n_loops, 1, even_edges) # output is always type 1
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        if not self.is_valid():
            return 0
        target_dim_sort = self.target.get_sort_dim()
        if target_dim_sort == 0:
            return 0
        return self.domain.n_edges * math.log(target_dim_sort, 2)

    def get_type(self):
        return 'contract edges'

    def operate_on(self, G):
        return self.ocontract.operate_on(G)


