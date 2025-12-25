"""Graph complexes based on ordinary (simple, without multiple edges) graphs.
Implemented Differentials: Contract edges, delete edges.
The generators only produce 1-vertex irreducible graphs."""


__all__ = ['graph_type', 'sub_types', 'OrdinaryGVS', 'OrdinaryGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
           'DeleteEdgesGO', 'DeleteEdgesD', 'OrdinaryGC']

import copy
import itertools
from sage.all import *
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import NautyInterface
import Parameters
import GCDimensions
from itertools import product

graph_type = "ordinary"

sub_types = {True: "even_edges", False: "odd_edges"}


# ------- Graph Vector Space --------
class OrdinaryGVS(GraphVectorSpace.GraphVectorSpace):
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
        super().__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops

    def __hash__(self):
        return hash("gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return (3 * self.n_vertices <= 2 * self.n_edges) and self.n_vertices > 0 and self.n_loops >= 0 \
            and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2

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
        return NautyInterface.list_simple_graphs(self.n_vertices, self.n_edges)

    def perm_sign(self, G, p):
        if self.even_edges:
            # The sign is (induced sign on vertices) * (induced sign edge orientations)
            sign = Shared.Perm(p).signature()
            for (u, v) in G.edges(labels=False,sort=True):
                # We assume the edge is always directed from the larger to smaller index
                if (u < v and p[u] > p[v]) or (u > v and p[u] < p[v]):
                    sign *= -1
            return sign
        else:
            # The sign is (induced sign of the edge permutation)
            # We assume the edges are always lexicographically ordered
            # For the computation we use that G.edges() returns the edges in lex ordering
            # We first label the edges on a copy of G lexicographically
            G1 = copy(G)
            Shared.enumerate_edges(G1)
            # We permute the graph, and read of the new labels
            G1.relabel(p, inplace=True)
            return Shared.Perm([j for (u, v, j) in G1.edges(sort=True)]).signature()

class LieBracket():
    @staticmethod
    def insertion_product(G1, G2, v, even_edges): 
        """Constructs the linear combination of graphs obtained by inserting G2 into vertex v of G1

        :param G1: The first graph
        :type G1: graph
        :param G2: The second graph
        :type G2: graph
        :param v: The vertex in G1 where G2 is inserted
        :type v: int
        """
        # sanity checks
        V1 = OrdinaryGVS(G1.order(), G1.size() - G1.order() + 1, even_edges)
        
        # permute the vertices of G1 so that v becomes the last vertex
        p1 = list(range(G1.order()))
        p1[v], p1[-1] = p1[-1], p1[v]
        G1b = copy(G1)
        G1b.relabel(p1, inplace=True)
        sgn = V1.perm_sign(G1, p1)

        # prepare G2 by shifting its vertex labels
        shift = G1.order() - 1
        p2 = [i + shift for i in range(G2.order())]
        G2b = copy(G2)
        G2b.relabel(p2, inplace=True)

        # remember neighbors of v in G1
        neighbors = list(G1b.neighbors(G1b.order() - 1))
        # remove v from G1b and take union with G2b
        G1b.delete_vertex(G1b.order() - 1)
        G = G1b.union(G2b)

        # now we have to reconnect the edges that were adjacent to v
        # we have to consider all possible ways to connect these edges to vertices of G2

        images = []
        for reconnection in itertools.product(range(G2.order()), repeat=len(neighbors)):
            Gc = copy(G)
            for i, u in enumerate(neighbors):
                v2 = reconnection[i] + shift
                Gc.add_edge(u, v2)
            
            images.append((Gc, sgn ))
        return images
    
    @staticmethod
    def lie_bracket_single(G1, G2, even_edges):
        """Computes the Lie bracket [G1, G2] where G1 is in V1 and G2 is in V2

        :param V1: The vector space of G1
        :type V1: OrdinaryGVS
        :param V2: The vector space of G2
        :type V2: OrdinaryGVS
        :param G1: The first graph
        :type G1: graph
        :param G2: The second graph
        :type G2: graph
        """
        V1 = OrdinaryGVS(G1.order(), G1.size() - G1.order() + 1, even_edges)
        V2 = OrdinaryGVS(G2.order(), G2.size() - G2.order() + 1, even_edges)
        
         # first term: insert G2 into G1
        images = []
        for v in range(G1.order()):
            images += OrdinaryGVS.insertion_product(G1, G2, v, even_edges)
        
        # relative sign of second term
        gsgn = -1
        if V1.even_edges:
            if (G1.size() * G2.size()) % 2 == 1:
                gsgn = 1
        else:
            if ( (G1.order()+1) * (G2.order()+1)) % 2 == 1:
                gsgn = 1
        

        for v in range(G2.order()):
            images_v2 = V2.insertion_product(G2, G1, v)
            for (G, sgn) in images_v2:
                images.append((G, gsgn * sgn))  # minus sign for second term in Lie bracket
        
        return images
    
    @staticmethod
    def lie_bracket_single_vector(G1, G2, even_edges):
        """Computes the Lie bracket [G1, G2] where G1 is in V1 and G2 is in V2

        :param V1: The vector space of G1
        :type V1: OrdinaryGVS
        :param V2: The vector space of G2
        :type V2: OrdinaryGVS
        :param G1: The first graph
        :type G1: graph
        :param G2: The second graph
        :type G2: graph
        """
        images = OrdinaryGVS.lie_bracket_single(G1, G2, even_edges)
        # now express in basis
        V1 = OrdinaryGVS(G1.order(), G1.size() - G1.order() + 1, even_edges)
        V2 = OrdinaryGVS(G2.order(), G2.size() - G2.order() + 1, even_edges)
        V = OrdinaryGVS(G1.order() + G2.order() - 1, G1.size() + G2.size() - (G1.order() + G2.order()) + 1, even_edges)
        basis_dict = V.get_g6_coordinates_dict()
        result = [0 for _ in range(len(basis_dict))]
        for (G, sgn) in images:
            (g6, sgn2) = V.graph_to_canon_g6(G)
            sgn_total = sgn * sgn2
            if g6 in basis_dict:
                index = basis_dict[g6]
                result[index] += sgn_total
        return result

    @staticmethod
    def lie_bracket_multi_vector(G1_list, G2_list, even_edges):
        """Computes the Lie bracket 
        """
        # first compute all pairwise Lie brackets
        images = []
        G1_first , _ = G1_list[0]
        G2_first , _ = G2_list[0]
        V1 = OrdinaryGVS(G1_first.order(), G1_first.size() - G1_first.order() + 1, even_edges)
        V2 = OrdinaryGVS(G2_first.order(), G2_first.size() - G2_first.order() + 1, even_edges)
        V = OrdinaryGVS(G1_first.order() + G2_first.order() - 1, G1_first.size() + G2_first.size() - (G1_first.order() + G2_first.order()) + 1, even_edges)
        
        basis_dict = V.get_g6_coordinates_dict()
        result = [0 for _ in range(len(basis_dict))]
        
        for (G1,x1) in G1_list:
            for (G2,x2) in G2_list:
                for G,y in OrdinaryGVS.lie_bracket_single(G1, G2, even_edges):
                    sgn_total = x1 * x2 * y
                    (g6, sgn2) = V.graph_to_canon_g6(G)
                    sgn_total *= sgn2
                    if g6 in basis_dict:
                        index = basis_dict[g6]
                        result[index] += sgn_total
        return result
        
        


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
        super().__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s_%s" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


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
        domain = OrdinaryGVS(n_vertices, n_loops, even_edges)
        target = OrdinaryGVS(n_vertices - 1, n_loops, even_edges)
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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image = []
        for (i, e) in enumerate(G.edges(labels=False,sort=True)):
            (u, v) = e
            # print("contract", u, v)
            pp = Shared.permute_to_left(
                (u, v), range(self.domain.n_vertices))
            sgn = self.domain.perm_sign(G, pp)
            G1 = copy(G)
            G1.relabel(pp, inplace=True)
            Shared.enumerate_edges(G1)
            previous_size = G1.size()
            G1.merge_vertices([0, 1])
            if (previous_size - G1.size()) != 1:
                continue
            # print(sgn)
            G1.relabel(list(range(G1.order())), inplace=True)
            if not self.domain.even_edges:
                # p = [j for (a, b, j) in G1.edges()]
                # sgn *= Permutation(p).signature()
                sgn *= Shared.shifted_edge_perm_sign(G1)
            else:
                sgn *= -1  # TODO overall sign for even edges
            image.append((G1, sgn))
        return image


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: OrdinaryGraphSumVS
        """
        super().__init__(sum_vector_space,
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
        super().__init__(domain, target)

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
        super().__init__(sum_vector_space,
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
class OrdinaryGC(GraphComplex.GraphComplex):
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

        sum_vector_space = OrdinaryGraphSumVS(v_range, l_range, even_edges,
                                              shift_loops_minus_vertices=shift_loops_minus_vertices)
        differential_list = []
        if not set(differentials) <= {'contract', 'delete'}:
            raise ValueError(
                "Differentials for ordinary graph complex: 'contract', 'delete'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)
        if 'delete' in differentials:
            delete_edges_dif = DeleteEdgesD(sum_vector_space)
            differential_list.append(delete_edges_dif)
        super().__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))
