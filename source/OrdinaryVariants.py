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
import OrdinaryGraphComplex


graph_type = "ordinary_variants"

sub_types = {True: "even_edges", False: "odd_edges"}


# ------- Graph Vector Space --------
class OrdinaryGVSFull(GraphVectorSpace.GraphVectorSpace):
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

    def __init__(self, n_vertices, n_loops, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, even_edges)
        super(OrdinaryGVSFull, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops

    def __hash__(self):
        return hash("gra_full%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra_full%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra_full%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return self.ogvs.is_valid()

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        return self.ogvs.get_work_estimate()

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        return NautyInterface.list_simple_graphs(self.n_vertices, self.n_edges, onlyonevi=False)

    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G, p)
    
# ------- Graph Vector Space --------
class OrdinaryGVSBridgeless(GraphVectorSpace.GraphVectorSpace):
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

    def __init__(self, n_vertices, n_loops, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, even_edges)
        self.full_ogvs = OrdinaryGVSFull(n_vertices, n_loops, even_edges)
        super(OrdinaryGVSBridgeless, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops

    def __hash__(self):
        return hash("gra_bl%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra_bl%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra_bl%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return self.ogvs.is_valid()

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        return self.ogvs.get_work_estimate()

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        for G in self.full_ogvs.get_basis():
            # check if sage graph is bridgeless
            if len(list(G.bridges())) == 0:
                yield G

    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G, p)
    

class OrdinaryGVSPanzer(GraphVectorSpace.GraphVectorSpace):
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

    def __init__(self, n_vertices, n_loops, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, even_edges)
        self.full_ogvs = OrdinaryGVSFull(n_vertices, n_loops, even_edges)
        super(OrdinaryGVSPanzer, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops

    def __hash__(self):
        return hash("gra_panzer%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra_panzer%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra_panzer%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return self.ogvs.is_valid()

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        return self.ogvs.get_work_estimate()

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        for G in self.full_ogvs.get_basis():
            if not G.is_triconnected():
                continue
            # check if the graph has a subgraph of negative degre
            #if 2*G.order() - G.size() - 2 < 0:
            #    continue
            ok = True
            for H in G.connected_subgraph_iterator():
                # only consider proper subgraphs
                if H.order() == G.order():
                    continue
                if 2*H.order() - H.size() - 2 < 0:
                    ok = False
                    break
            if ok:
                yield G

    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G, p)
    




class OrdinaryGVSTriconnected(GraphVectorSpace.GraphVectorSpace):
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

    def __init__(self, n_vertices, n_loops, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, even_edges)
        self.full_ogvs = OrdinaryGVSFull(n_vertices, n_loops, even_edges)
        super(OrdinaryGVSTriconnected, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops

    def __hash__(self):
        return hash("gra_tri%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra_tri%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra_tri%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return self.ogvs.is_valid()

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        return self.ogvs.get_work_estimate()

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        for G in self.ogvs.get_basis():
            # check if sage graph is triconnected
            if G.is_triconnected():
                yield G

    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G, p)
    
    
class OrdinaryGVSKVconnected(GraphVectorSpace.GraphVectorSpace):
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

    def __init__(self, n_vertices, n_loops, k_conn, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.k_conn = k_conn
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, even_edges)
        self.full_ogvs = OrdinaryGVSFull(n_vertices, n_loops, even_edges)
        super(OrdinaryGVSKVconnected, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.k_conn == other.k_conn

    def __hash__(self):
        return hash("gra_kv%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra_kv%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    # def get_ref_basis_file_path(self):
    #     s = "gra_tri%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('vertconnectivity',self.k_conn)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return self.ogvs.is_valid()

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        return self.ogvs.get_work_estimate()

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        for G in self.full_ogvs.get_basis():
            # check if sage graph is triconnected
            if G.vertex_connectivity() >= self.k_conn:
                yield G

    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G, p)
    

class OrdinaryGVSKEconnected(GraphVectorSpace.GraphVectorSpace):
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

    def __init__(self, n_vertices, n_loops, k_conn, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.k_conn = k_conn
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(n_vertices, n_loops, even_edges)
        self.full_ogvs = OrdinaryGVSFull(n_vertices, n_loops, even_edges)
        super(OrdinaryGVSKEconnected, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.k_conn == other.k_conn

    def __hash__(self):
        return hash("gra_ke%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra_ke%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    # def get_ref_basis_file_path(self):
    #     s = "gra_tri%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('edgeconnectivity',self.k_conn)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return self.ogvs.is_valid()

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        return self.ogvs.get_work_estimate()

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        for G in self.full_ogvs.get_basis():
            # check if sage graph is triconnected
            if G.edge_connectivity() >= self.k_conn:
                yield G

    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G, p)


class OrdinaryGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of ordinary graph vector spaces with specified edge parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, even_edges, shift_loops_minus_vertices=1):
        """Initialize the sum vector space.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.v_range = v_range
        self.l_range = l_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        if shift_loops_minus_vertices is None:
            vs_list = [OrdinaryGVS(v, l, self.even_edges) for (
                v, l) in itertools.product(self.v_range, self.l_range)]
        else:
            vs_list = []
            for (v, l) in itertools.product(self.v_range, self.l_range):
                if l - v <= shift_loops_minus_vertices:
                    vs_list.append(OrdinaryGVS(v, l, self.even_edges))
        super(OrdinaryGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s_%s" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


# ------- Operators --------
class ContractEdgesGOFull(GraphOperator.GraphOperator):
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
        if not ContractEdgesGOFull.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        self.ocontract = OrdinaryGraphComplex.ContractEdgesGO(domain.ogvs, target.ogvs)
        super(ContractEdgesGOFull, self).__init__(domain, target)

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
    def generate_operator(cls, n_vertices, n_loops, even_edges):
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
        domain = OrdinaryGVSFull(n_vertices, n_loops, even_edges)
        target = OrdinaryGVSFull(n_vertices - 1, n_loops, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD_full%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD_full%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD_full%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD_full%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        return self.ocontract.operate_on(G)


class ContractEdgesGOBridgeless(GraphOperator.GraphOperator):
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
        if not ContractEdgesGOBridgeless.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        self.ocontract = OrdinaryGraphComplex.ContractEdgesGO(domain.ogvs, target.ogvs)
        super(ContractEdgesGOBridgeless, self).__init__(domain, target)

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
    def generate_operator(cls, n_vertices, n_loops, even_edges):
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
        domain = OrdinaryGVSBridgeless(n_vertices, n_loops, even_edges)
        target = OrdinaryGVSBridgeless(n_vertices - 1, n_loops, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD_bl%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD_bl%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD_bl%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD_bl%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        return self.ocontract.operate_on(G)


class ContractEdgesGOPanzer(GraphOperator.GraphOperator):
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
        if not ContractEdgesGOPanzer.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        self.ocontract = OrdinaryGraphComplex.ContractEdgesGO(domain.ogvs, target.ogvs)
        super(ContractEdgesGOPanzer, self).__init__(domain, target)

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
    def generate_operator(cls, n_vertices, n_loops, even_edges):
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
        domain = OrdinaryGVSPanzer(n_vertices, n_loops, even_edges)
        target = OrdinaryGVSPanzer(n_vertices - 1, n_loops, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD_panzer%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD_panzer%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD_panzer%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD_panzer%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        return self.ocontract.operate_on(G)



class ContractEdgesGOTriconnected(GraphOperator.GraphOperator):
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
        if not ContractEdgesGOTriconnected.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        self.ocontract = OrdinaryGraphComplex.ContractEdgesGO(domain.ogvs, target.ogvs)
        super(ContractEdgesGOTriconnected, self).__init__(domain, target)

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
    def generate_operator(cls, n_vertices, n_loops, even_edges):
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
        domain = OrdinaryGVSTriconnected(n_vertices, n_loops, even_edges)
        target = OrdinaryGVSTriconnected(n_vertices - 1, n_loops, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD_tri%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD_tri%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD_tri%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD_tri%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        return self.ocontract.operate_on(G)


class ContractEdgesGOKV(GraphOperator.GraphOperator):
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
        if not ContractEdgesGOKV.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        self.ocontract = OrdinaryGraphComplex.ContractEdgesGO(domain.ogvs, target.ogvs)
        super(ContractEdgesGOKV, self).__init__(domain, target)

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
            and domain.even_edges == target.even_edges and domain.k_conn == target.k_conn

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, k_conn, even_edges):
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
        domain = OrdinaryGVSKVconnected(n_vertices, n_loops, k_conn, even_edges)
        target = OrdinaryGVSKVconnected(n_vertices - 1, n_loops, k_conn, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD_kv%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD_kv%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    # def get_ref_matrix_file_path(self):
    #     s = "contractD_tri%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    # def get_ref_rank_file_path(self):
    #     s = "contractD_tri%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        return self.ocontract.operate_on(G)

class ContractEdgesGOKE(GraphOperator.GraphOperator):
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
        if not ContractEdgesGOKE.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        self.ocontract = OrdinaryGraphComplex.ContractEdgesGO(domain.ogvs, target.ogvs)
        super(ContractEdgesGOKE, self).__init__(domain, target)

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
            and domain.even_edges == target.even_edges and domain.k_conn == target.k_conn

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, k_conn, even_edges):
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
        domain = OrdinaryGVSKEconnected(n_vertices, n_loops, k_conn, even_edges)
        target = OrdinaryGVSKEconnected(n_vertices - 1, n_loops, k_conn, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD_ke%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD_ke%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    # def get_ref_matrix_file_path(self):
    #     s = "contractD_tri%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    # def get_ref_rank_file_path(self):
    #     s = "contractD_tri%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        return self.ocontract.operate_on(G)

class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: OrdinaryGraphSumVS
        """
        super(ContractEdgesD, self).__init__(sum_vector_space,
                                             ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class DeleteEdgesGO(GraphOperator.GraphOperator):
    """Delete edges graph operator.

    Operates on an ordinary graph by deleting an edge.
    Only for graphs with odd edges.

    Attributes:
        - sub_type (str): Graphs sub type of the domain.
    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the delete edges graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: OrdinaryGVS
        :param target: Target vector space of the operator.
        :type target: OrdinaryGVS
        """
        if not DeleteEdgesGO.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for delete edges operator")
        self.sub_type = domain.sub_type
        super(DeleteEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding delete edges graph operator.

        The delete edge operator reduces the number of loops by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: OrdinaryGVS
        :param target: Potential target vector space of the operator.
        :type target: OrdinaryGVS
        :return: True if domain and target match to generate a corresponding contract edges graph operator.
        :rtype: bool
        """
        return domain.n_vertices == target.n_vertices and domain.n_loops - 1 == target.n_loops \
            and domain.even_edges == target.even_edges

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, even_edges):
        """Return a delete edge graph operator.

        :param n_vertices: Number of vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :return: Contract edges graph operator based on the specified domain vector space.
        :rtype: ContractEdgesGO
        """
        domain = OrdinaryGVS(n_vertices, n_loops, even_edges)
        target = OrdinaryGVS(n_vertices, n_loops - 1, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "delete_edgesD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "delete_edgesD%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "delete_edgesD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "delete_edgesD%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        return 'delete edges'

    def operate_on(self, G):
        # Operates on the graph G by deleting an edge.
        image = []
        for (i, e) in enumerate(G.edges(labels=False,sort=True)):
            (u, v) = e
            G1 = copy(G)
            G1.delete_edge((u, v))
            sgn = -1 if i % 2 else 1
            image.append((G1, sgn))
        return image


class DeleteEdgesD(GraphOperator.Differential):
    """Delete edges differential.

    Only for graphs with odd edges.
    """

    def __init__(self, sum_vector_space):
        """Initialize the delete edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: OrdinaryGraphSumVS
        """
        super(DeleteEdgesD, self).__init__(sum_vector_space,
                                           DeleteEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'delete edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_delete_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_delete_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_delete_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


# ------- Graph Complexes --------
class OrdinaryGCFull(GraphComplex.GraphComplex):
    """Graph complex for ordinary graphs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, even_edges, differentials, shift_loops_minus_vertices=1):
        """Initialize the graph complex.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param differentials: List of differentials. Options: 'contract', 'delete'.
        :type differentials: list(str)
        """
        self.v_range = v_range
        self.l_range = l_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        sum_vector_space = OrdinaryGraphSumVSFull(v_range, l_range, even_edges,
                                              shift_loops_minus_vertices=shift_loops_minus_vertices)
        differential_list = []
        if not set(differentials) <= {'contract', 'delete'}:
            raise ValueError(
                "Differentials for ordinary graph complex: 'contract', 'delete'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesDFull(sum_vector_space)
            differential_list.append(contract_edges_dif)
        if 'delete' in differentials:
            delete_edges_dif = DeleteEdgesD(sum_vector_space)
            differential_list.append(delete_edges_dif)
        super(OrdinaryGCFull, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))
