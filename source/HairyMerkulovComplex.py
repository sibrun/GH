"""Hairy graph complex, with vertices of valence 3 and 4.
Merkulov complex. """


__all__ = ['graph_type', 'sub_types', 'OrdinaryGVS', 'OrdinaryGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
           'DeleteEdgesGO', 'OrdinaryGC']

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
import os


graph_type = "hairyme"

sub_types = {(True, True): "even_edges_even_hairs", (True, False): "even_edges_odd_hairs",
             (False, True): "odd_edges_even_hairs", (False, False): "odd_edges_odd_hairs"}

# Option to include zero hairs in the hairy graph complexes.
zero_hairs = False


# ------- Graph Vector Space --------
class HairyMerkulovGVS(GraphVectorSpace.GraphVectorSpace):
    """Hairy graph vector space.
    Vertices have valence 3-6, with at most one vertex of valence 5 or 6

    Attributes:
        - n_vertices (int): Number of vertices.
        - n_loops (int): Number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - n_edges (int): Number of edges.
        - sub_type (str): Sub type of graphs.

    """

    def __init__(self, n_vertices, n_loops, n_hairs, even_edges, even_hairs, valence_type):
        """Initialize the ordinary graph vector space with (almost) all vertices of valence 3 and 4.

        :param n_vertices: int: Number of vertices.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param valence_type: Can be 34 (only 3 and 4-valent vertices),
            3456 (at most one vertex of valence 5 or 6), 56 (exactly one vertex of valence 5 or 6).
        """
        self.n_vertices = n_vertices
        self.n_hairs = n_hairs
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get( (self.even_edges, self.even_hairs) )
        self.valence_type = valence_type
        if valence_type not in { 34, 3456, 56 }:
            raise ValueError(f"OrdinaryMerkulovGVS: valence type must be in 34, 3456, 56")
        self.ogvs : HairyGraphComplex.HairyGraphVS = HairyGraphComplex.HairyGraphVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs)
        if valence_type != 3456:
            self.gvs3456 = HairyMerkulovGVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs, 3456)
        super().__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.even_edges == other.even_edges and self.valence_type == other.valence_type

    def __hash__(self):
        return hash("megra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "megra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('valence_type', self.valence_type), ('vertices', self.n_vertices), ('loops', self.n_loops), ('hairs', self.n_hairs)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        if self.valence_type == 34 and self.n_edges > 2*self.n_vertices - self.n_hairs/2:
            return False
        if self.valence_type == 56 and 2*self.n_edges < 3*self.n_vertices + 2 - self.n_hairs:
            return False
        return (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs) \
            and self.n_vertices > 0 \
            and self.n_loops >= 0 \
            and ((self.n_hairs >= 0) if zero_hairs else (self.n_hairs > 0)) \
            and self.n_edges <= 2*self.n_vertices+1 - self.n_hairs \
            and self.n_vertices >= self.n_hairs \
            and self.n_edges <= self.n_vertices * (self.n_vertices-1) /2

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return GCDimensions.get_hairy_dim_estimate(self.n_vertices, self.n_loops, self.n_hairs)
        #return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return

        n_vertices_1 = self.n_vertices
        n_vertices_2 = self.n_hairs + self.n_edges
        n_edges_bip = self.n_hairs + 2 * self.n_edges
        deg_range_2 = (1, 2)

        if self.valence_type == 3456:
            deg_range_1 = (3, 6)
            bipartite_graphs = ( self.ogvs._bip_to_ordinary(G)
                    for G in NautyInterface.list_bipartite_graphs(
                    n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip) )
            for G in bipartite_graphs:
                # check there is at most one vertex of valence >4
                if sum( (1 if len(G[v])>4 else 0) for v in range(self.n_vertices) ) <= 1:
                    yield G
        elif self.valence_type == 34:
            deg_range_1 = (3, 4)
            bipartite_graphs = ( self.ogvs._bip_to_ordinary(G)
                    for G in NautyInterface.list_bipartite_graphs(
                    n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip) )
            for G in bipartite_graphs:
                yield G
        elif self.valence_type == 56:
            deg_range_1 = (3, 6)
            bipartite_graphs = ( self.ogvs._bip_to_ordinary(G)
                    for G in NautyInterface.list_bipartite_graphs(
                    n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip) )
            for G in bipartite_graphs:
                # check there is exactly one vertex of valence >4
                if sum( (1 if len(G[v])>4 else 0) for v in range(self.n_vertices) ) == 1:
                    yield G

    def get_partition(self):
        return self.ogvs.get_partition()

    def perm_sign(self, G, p):
        return self.ogvs.perm_sign(G,p)


class HairyMerkulovGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of ordinary graph vector spaces with specified edge parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_range, even_edges, even_hairs, valence_types):
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
        self.h_range = h_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = sub_types.get(self.even_edges)

        vs_list = [HairyMerkulovGVS(v, l, h, self.even_edges, self.even_hairs, vt)
                for v in self.v_range
                for l in self.l_range
                for h in self.h_range
                for vt in valence_types ]
        super().__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('hairs', self.h_range)])

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

    def __init__(self, domain:HairyMerkulovGVS, target:HairyMerkulovGVS):
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
        self.oop = HairyGraphComplex.ContractEdgesGO(domain.ogvs, target.ogvs)
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
        return domain.n_vertices - 1 == target.n_vertices \
            and domain.n_loops == target.n_loops \
            and domain.n_hairs == target.n_hairs \
            and domain.even_edges == target.even_edges \
            and domain.even_hairs == target.even_hairs \
            and domain.valence_type == 34 and target.valence_type in { 3456, 56 }

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, even_edges, even_hairs, to3456 = True):
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
        domain = HairyMerkulovGVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs, 34)
        target = HairyMerkulovGVS(n_vertices - 1, n_loops, n_hairs, even_edges, even_hairs, 3456 if to3456 else 56)

        return cls(domain, target)

    def get_matrix_file_path(self):
        s = f"contractD{self.target.valence_type}_{self.domain.n_vertices}_{self.domain.n_loops}_{self.domain.n_hairs}.txt"
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = f"contractD{self.target.valence_type}_{self.domain.n_vertices}_{self.domain.n_loops}_{self.domain.n_hairs}_rank.txt"
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

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
        lst = self.oop.operate_on(G)
        # g6, _ = self.domain.ogvs.graph_to_canon_g6(G)
        # print(g6," ->")
        # for (GG,v) in lst:
        #     print(GG.graph6_string())
        return lst


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


def cohom_formatted_merkulov(D1, D2, Dc2):
    vs = D1.get_domain()
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    rc2 = 0
    if D1.is_valid():
        if D1.exists_rank_file():
            r1 = D1.get_matrix_rank()
        else:
            return "?"
    if D2.is_valid():
        if D2.exists_rank_file():
            r2 = D2.get_matrix_rank()
        else:
            return "?"

    if Dc2.is_valid():
        if Dc2.exists_rank_file():
            rc2 = Dc2.get_matrix_rank()
        else:
            return "?"

    # exact or not?
    r_str = "" if D1.exists_exact_rank() and D2.exists_exact_rank(
    ) and Dc2.exists_exact_rank() else " p"

    return str(d+rc2-r1-r2) + r_str

def get_34cohom_dim(v, l, h, even_e, even_h):
    """ Compute cohomology dimension ..."""
    op1 = ContractEdgesGO.generate_operator(v, l, h, even_e, even_h)
    op2 = ContractEdgesGO.generate_operator(v + 1, l, h, even_e, even_h)
    opc = ContractEdgesGO.generate_operator(v + 1, l, h, even_e, even_h, False)

    return cohom_formatted_merkulov(op1, op2, opc)

    # return vs34.get_34dimension() - D34rank -DD34rank + DD5rank

def cohom_formatted2(D1, D2):
    vs = D1.get_domain()
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    if D1.is_valid():
        if D1.exists_rank_file():
            r1 = D1.get_matrix_rank()
        else:
            return "?"
    if D2.is_valid():
        if D2.exists_rank_file():
            r2 = D2.get_matrix_rank()
        else:
            return "?"

    # exact or not?
    r_str = "" if D1.exists_exact_rank() and D2.exists_exact_rank() else " p"

    return str(d-r1-r2) + r_str

def get_ref_cohom_dim(v, l, h, even_e, even_h):
    """ Compute cohomology dimension ..."""
    op1 = HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, even_e, even_h)
    op2 = HairyGraphComplex.ContractEdgesGO.generate_operator(v+1, l, h, even_e, even_h)

    return cohom_formatted2(op1, op2)

# ------- Graph Complexes --------
class HairyMerkulovGC(GraphComplex.GraphComplex):
    """Graph complex for ordinary graphs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_range, even_edges, even_hairs, differentials):
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
        self.h_range = h_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = sub_types.get(self.even_edges)

        sum_vector_space = HairyMerkulovGraphSumVS(v_range, l_range, h_range, even_edges, even_hairs, [34, 3456, 56] )
        differential_list = []
        if not set(differentials) <= {'contract'}:
            raise ValueError(
                "Differentials for ordinary graph complex: 'contract'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)

        super().__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))

    def print_cohom(self):
        for h in self.h_range:
            for l in self.l_range:
                print(f"l={l}, h={h}: {[ get_34cohom_dim(v,l,h, self.even_edges, self.even_hairs) for v in self.v_range]}")


    def print_cohom_reference(self):
        """Print cohomology of analogous ordinary hairy gvs for reference."""
        for h in self.h_range:
            for l in self.l_range:
                print(f"l={l}, h={h}: {[ get_ref_cohom_dim(v,l,h, self.even_edges, self.even_hairs) for v in self.v_range]}")


