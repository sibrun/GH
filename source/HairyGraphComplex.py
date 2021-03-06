"""Graph complexes based on simple graphs with hairs. Without multiple edges and multiple hairs per vertex.
One hair is composed of a hair vertex and an edge connecting it to a vertex. The parity of the hair refers to the parity
of the hair vertex alone.
Implemented Differentials: Contract edges, edge to one hair."""


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

graph_type = "hairy"

sub_types = {(True, True): "even_edges_even_hairs", (True, False): "even_edges_odd_hairs",
             (False, True): "odd_edges_even_hairs", (False, False): "odd_edges_odd_hairs"}

zero_hairs = False      # Option to include zero hairs in the hairy graph complexes.


# ------- Graph Vector Space --------
class HairyGraphVS(GraphVectorSpace.GraphVectorSpace):
    """Hairy graph vector space.

    Sub vector space with specified number of vertices, loops, hairs, even or odd edges, even or odd hair vertices
    and at least trivalent vertices. No multiple edges and not mor than one hair is attached to a vertex. One hair is
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

    def __init__(self, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        """Initialize the hairy graph vector space.

        :param n_vertices: Number of internal vertices.
        :type n_vertices: int
        :param n_loops: Number of loops.
        :type n_loops: int
        :param n_hairs: Number of hairs.
        :type n_hairs: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs: Parity of the hair vertices. True for even hairs, False for odd hairs.
        :type even_hairs: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs = n_hairs
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))
        super(HairyGraphVS, self).__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(self.n_vertices + self.n_hairs, self.n_loops, self.even_edges)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs

    def get_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('hairs', self.n_hairs)])

    def get_partition(self):
        # All internal vertices are in color 1, the hair vertices are in color 2.
        return [list(range(0, self.n_vertices)), list(range(self.n_vertices, self.n_vertices + self.n_hairs))]

    def is_valid(self):
        # At least trivalent internal vertices.
        l = (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # Positive number of vertices, non negative number of loops, non-negative or positive number of hairs.
        l = l and self.n_vertices > 0 and self.n_loops >= 0 and \
            ((self.n_hairs >= 0) if zero_hairs else (self.n_hairs > 0))
        # At most a full graph.
        l = l and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2
        # At most one hair per vertex.
        l = l and self.n_vertices >= self.n_hairs
        return l

    def get_work_estimate(self):
        # TODO
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_generating_graphs(self):
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2.
        # Degree 1 vertices are hairs, degree 2 vertices are edges and are removed later.
        # No multiple hairs and edges.
        if not self.is_valid():
            return []
        n_vertices_1 = self.n_vertices
        n_vertices_2 = self.n_hairs + self.n_edges
        n_edges_bip = self.n_hairs + 2 * self.n_edges
        deg_range_1 = (3, n_edges_bip + 1)
        deg_range_2 = (1, 2)
        bipartite_graphs = NautyInterface.list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip)
        return [self._bip_to_ordinary(G) for G in bipartite_graphs]

    def perm_sign(self, G, p):
        # The sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the hair-vertices.
        sgn = self.ogvs.perm_sign(G, p)
        # Compute the extra contribution from hairs.
        if self.even_hairs == self.even_edges:
            hairs = p[self.n_vertices:]
            if len(hairs) != 0:
                sgn *= Shared.Perm.shifted(hairs).signature()
        return sgn

    def _bip_to_ordinary(self, G):
        # Translates bipartite into ordinary graph by replacing a bivalent vertex of colour 2 with an edge.
        for v in range(self.n_vertices, self.n_vertices + self.n_hairs + self.n_edges):
            neighbors = G.neighbors(v)
            n_l = len(neighbors)
            if n_l == 1: #hair
                continue
            elif n_l == 2: #edge
                G.add_edge(neighbors)
                G.delete_vertex(v)
            else:
                raise ValueError('%s: Vertices of second colour should have 1 or 2 neighbours' % str(self))
        return G


class HairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of hairy graph vector spaces with specified edge and hair parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_range (range): Range for the number of hairs.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs (bool): True for even hairs, False for odd hairs.
        - sub_type (str): Sub type of graphs.
    """
    def __init__(self, v_range, l_range, h_range, even_edges, even_hairs):
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
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))

        vs_list = [HairyGraphVS(v, l, h, self.even_edges, self.even_hairs) for
                   (v, l, h) in itertools.product(self.v_range, self.l_range, self.h_range)]
        super(HairyGraphSumVS, self).__init__(vs_list)

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

    Operate on a hairy graph by contracting an edge not connected to a hair vertex and unifying the two adjacent vertices.

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
        super(ContractEdgesGO, self).__init__(domain, target)

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
               and domain.n_hairs == target.n_hairs and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
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
        domain = HairyGraphVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs)
        target = HairyGraphVS(n_vertices - 1, n_loops, n_hairs, even_edges, even_hairs)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        if not self.is_valid():
            return 0
        try:
            (domain_dim, dimtarget_dim) = (self.domain.get_dimension(), self.target.get_dimension())
        except StoreLoad.FileNotFoundError:
            return 0
        if domain_dim == 0 or dimtarget_dim == 0:
            return 0
        return self.domain.n_edges * math.log(self.target.get_dimension(), 2)

    def get_type(self):
        return 'contract edges'

    def operate_on(self,G):
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # only edges not connected to a hair-vertex can be contracted
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices:
                continue
            pp = Shared.permute_to_left((u, v), range(0, self.domain.n_vertices + self.domain.n_hairs))
            sgn = self.domain.perm_sign(G, pp)
            G1 = copy(G)
            G1.relabel(pp, inplace=True)
            Shared.enumerate_edges(G1)
            previous_size = G1.size()
            G1.merge_vertices([0, 1])
            if (previous_size - G1.size()) != 1:
                continue
            G1.relabel(list(range(0, G1.order())), inplace=True)
            if not self.domain.even_edges:
                sgn *= Shared.shifted_edge_perm_sign(G1)
            image.append((G1, sgn))
        return image


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""
    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: HairyGraphSumVS
        """
        super(ContractEdgesD, self).__init__(sum_vector_space, ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class EdgeToOneHairGO(GraphOperator.GraphOperator):
    """Edge to one hair graph operator.

    Operates on a hairy graph by deleting an edge and adding a hair to one of the vertices adjacent to the
    deleted edge.
    Only for graphs with odd hairs.

    Attributes:
        - sub_type (str): Graphs sub type of the domain.
    """
    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the edge to one hair graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: HairyGraphVS
        :param target: Target vector space of the operator.
        :type target: HairyGraphVS
        """
        self.sub_type = domain.sub_type
        super(EdgeToOneHairGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding edge to one hair graph operator.

        The edge to one hair operator reduces the number of loops by one and increases the number of hairs by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: HairyGraphVS
        :param target: Potential target vector space of the operator.
        :type target: HairyGraphVS
        :return: True if domain and target match to generate a corresponding edge to one hair graph operator.
        :rtype: bool
        """
        return domain.n_vertices == target.n_vertices and domain.n_loops - 1 == target.n_loops \
               and domain.n_hairs + 1 == target.n_hairs and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        """Return an edge to one hair graph operator.

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
        :return: Edge to one hair graph operator based on the specified domain vector space.
        :rtype: EdgeToOneHairGO
        """
        domain = HairyGraphVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs)
        target = HairyGraphVS(n_vertices, n_loops - 1, n_hairs + 1, even_edges, even_hairs)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "edge_to_one_hairD%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "edge_to_one_hairD%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        if not self.is_valid():
            return 0
        try:
            target_dim = self.target.get_dimension()
        except StoreLoad.FileNotFoundError:
            return 0
        if target_dim == 0:
            return 0
        return self.domain.n_edges * math.log(self.target.get_dimension(), 2)

    def get_type(self):
        return 'edge to one hair'

    def operate_on(self,G):
        # Operate on a hairy graph by deleting an edge and adding a hair to one of the vertices adjacent to the
        # deleted edge.
        sgn0 = -1 if G.order() % 2 else 1
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # Only edges not connected to a hair-vertex can be cut
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices:
                continue
            G1 = copy(G)
            if not self.domain.even_edges:
                Shared.enumerate_edges(G1)
                e_label = G1.edge_label(u, v)
            G1.delete_edge((u, v))
            new_hair_idx = self.domain.n_vertices + self.domain.n_hairs
            G1.add_vertex(new_hair_idx)
            G2 = copy(G1)
            G1.add_edge((u, new_hair_idx))
            G2.add_edge((v, new_hair_idx))
            if not self.domain.even_edges:
                G1.set_edge_label(u, new_hair_idx, e_label)
                G2.set_edge_label(v, new_hair_idx, e_label)
                sgn1 = Shared.edge_perm_sign(G1)
                sgn2 = Shared.edge_perm_sign(G2)
            else:
                sgn1 = 1
                sgn2 = -1
            image.append((G1, sgn1*sgn0))
            image.append((G2, sgn2*sgn0))
        return image


class EdgeToOneHairD(GraphOperator.Differential):
    """Edge to one hair differential.

    Only for graphs with odd hairs.
    """
    def __init__(self, sum_vector_space):
        """Initialize the edge to one hair differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: HairyGraphSumVS
        """
        super(EdgeToOneHairD, self).__init__(sum_vector_space, EdgeToOneHairGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'edge to one hair'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_edge_to_one_hair_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_edge_to_one_hair_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_plot_parameter_order(self):
        return (1, 2, 0)


# ------- Graph Complex --------
class HairyGC(GraphComplex.GraphComplex):
    """Graph complex for hairy graphs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_range (range): Range for the number of hairs.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs (bool): True for even hairs, False for odd hairs.
        - sub_type (str): Sub type of graphs.
    """
    def __init__(self, v_range, l_range, h_range, even_edges, even_hairs, differentials):
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
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))

        sum_vector_space = HairyGraphSumVS(self.v_range, self.l_range, self.h_range, self.even_edges, self.even_hairs)
        differential_list = []
        if not differentials <= {'contract', 'et1h'}:
            raise ValueError("Differentials for hairy graph complex: 'contract', 'et1h'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)
        if 'et1h' in differentials:
            edge_to_one_hair_dif = EdgeToOneHairD(sum_vector_space)
            differential_list.append(edge_to_one_hair_dif)
        super(HairyGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))