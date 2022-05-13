"""Graph complexes based on simple graphs with two colours of hairs. Without multiple edges and multiple hairs of the
same colour per vertex.
Implemented Differentials: Contract edges, split edges."""


__all__ = ['graph_type', 'get_sub_type', 'BiColoredHairyGraphVS', 'BiColoredHairyGraphSumVS', 'ContractEdgesGO',
           'ContractEdgesD', 'SplitEdgesGO', 'SplitEdgesD', 'BiColoredHairyGC']

import itertools
from sage.all import *
from sage.combinat.shuffle import ShuffleProduct
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import NautyInterface
import OrdinaryGraphComplex
import HairyGraphComplex
import StoreLoad
import Parameters

graph_type = "bi_colored_hairy"


def get_sub_type(even_edges, even_hairs_a, even_hairs_b):
    sub_type = ('even' if even_edges else 'odd') + '_edges'
    sub_type += '_' + ('even' if even_hairs_a else 'odd') + '_hairs_a'
    sub_type += '_' + ('even' if even_hairs_b else 'odd') + '_hairs_b'
    return sub_type


# Option to include zero hairs in the hairy graph complexes.
zero_hairs = False


# ------- Graph Vector Space --------
class BiColoredHairyGraphVS(GraphVectorSpace.GraphVectorSpace):
    """Hairy graph vector space with two colours of hairs.

    Sub vector space with specified number of vertices, loops, hairs per colour, even or odd edges, even or odd hair vertices
    per colour and at least trivalent vertices. No multiple edges and not mor than one hair per colour is attached to a vertex.
    One hair is composed of a hair vertex and an edge connecting it to a vertex. The parity of the hair refers to the
    parity of the hair vertex alone.

    Attributes:
        - n_vertices (int): Number of internal vertices.
        - n_loops (int): Number of loops.
        - n_hairs_a (int): Number of hairs of the first colour.
        - n_hairs_b (int): Number of hairs of the second colour.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs_a (bool): Parity of the hair vertices of the first colour. True for even hairs_a, False for odd hairs_a.
        - even_hairs_b (bool): Parity of the hair vertices of the second colour. True for even hairs_b, False for odd hairs_b.
        - n_edges (int): Number of edges.
        - sub_type (str): Sub type of graphs.
        - ogvs (OrdinaryGraphComplex.OrdinaryGVS): Ordinary graph vector space without hairs.

    """

    def __init__(self, n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b):
        """Initialize the bi colored hairy graph vector space.

        :param n_vertices: Number of internal vertices.
        :type n_vertices: int
        :param n_loops: Number of loops.
        :type n_loops: int
        :param n_hairs_a: Number of hairs of the first colour.
        :type n_hairs_a: int
        :param n_hairs_b: Number of hairs of the second colour.
        :type n_hairs_b: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: Parity of the hair vertices of the first colour. True for even hairs_a, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: Parity of the hair vertices of the second colour. True for even hairs_a, False for odd hairs_b.
        :type even_hairs_b: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs_a = n_hairs_a
        self.n_hairs_b = n_hairs_b
        self.n_hairs = self.n_hairs_a + self.n_hairs_b
        self.even_edges = even_edges
        self.even_hairs_a = even_hairs_a
        self.even_hairs_b = even_hairs_b
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = get_sub_type(
            self.even_edges, self.even_hairs_a, self.even_hairs_b)
        super(BiColoredHairyGraphVS, self).__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(
            self.n_vertices + self.n_hairs, self.n_loops, self.even_edges)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and \
            self.n_hairs_a == other.n_hairs_a and self.n_hairs_b == other.n_hairs_b

    def __str__(self):
        return f"{graph_type} ({self.sub_type}) {self.n_vertices} vertices, {self.n_loops} loops, {self.n_hairs_a} hairs a, {self.n_hairs_b} hairs b"

    def __hash__(self):
        return hash(str(self))

    def get_basis_file_path(self):
        s = "gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops),
                                   ('hairs_a', self.n_hairs_a), ('hairs_b', self.n_hairs_b)])

    def get_partition(self):
        # All internal vertices are in color 1, the hairs_a vertices are in color 2, and the hairs_b vertices are in color 3.
        return [list(range(0, self.n_vertices)), list(range(self.n_vertices, self.n_vertices + self.n_hairs_a)),
                list(range(self.n_vertices + self.n_hairs_a, self.n_vertices + self.n_hairs))]

    def is_valid(self):
        # At least trivalent internal vertices
        l = (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # Positive number of vertices, non negative number of loops, non-negative or positive number of hairs for each color.
        l = l and self.n_vertices > 0 and self.n_loops >= 0 and \
            ((self.n_hairs_a >= 0 and self.n_hairs_b >= 0) if zero_hairs
             else (self.n_hairs_a > 0 and self.n_hairs_b > 0))
        # At most a full graph.
        l = l and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2
        # At most one hair of each color per vertex.
        l = l and self.n_vertices >= max(self.n_hairs_a, self.n_hairs_b)
        return l

    def get_work_estimate(self):
        # TODO
        # Return the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return (self.n_vertices ** self.n_hairs_a) * (self.n_vertices ** self.n_hairs_b) * binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / (factorial(self.n_vertices) * factorial(self.n_hairs_a) * factorial(self.n_hairs_b))

    def get_generating_graphs(self):
        # First produce all hairy graphs in the same way as for hairy graph vector spaces.
        # Then assign the hairs to colour a and b respectively by generating all shuffles of the hairs.
        if not self.is_valid():
            return []
        n_vertices_1 = self.n_vertices
        n_vertices_2 = self.n_hairs + self.n_edges
        n_edges_bip = self.n_hairs + 2 * self.n_edges
        deg_range_1 = (3, n_edges_bip)
        deg_range_2 = (1, 2)
        bipartite_graphs = NautyInterface.list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2,
                                                                n_edges_bip)

        return (GG for G in bipartite_graphs
                for GG in self._hairy_to_bi_colored_hairy(self._bip_to_ordinary(G))
                )
        # hairy_graphs = [self._bip_to_ordinary(G) for G in bipartite_graphs]
        # list_of_lists = [self._hairy_to_bi_colored_hairy(G) for G in hairy_graphs]
        # return list(itertools.chain.from_iterable(list_of_lists))

    def perm_sign(self, G, p):
        # The sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the hair-vertices.
        sgn = self.ogvs.perm_sign(G, p)
        # Compute the extra contribution from hairs_a and hairs_b separately.
        if self.even_hairs_a == self.even_edges:
            hairs_a = p[self.n_vertices: self.n_vertices + self.even_hairs_a]
            if len(hairs_a) != 0:
                sgn *= Shared.Perm.shifted(hairs_a).signature()
        if self.even_hairs_b == self.even_edges:
            hairs_b = p[self.n_vertices +
                        self.even_hairs_a: self.n_vertices + self.n_hairs]
            if len(hairs_b) != 0:
                sgn *= Shared.Perm.shifted(hairs_b).signature()
        return sgn

    def _bip_to_ordinary(self, G):
        # Translate bipartite into ordinary graph by replacing a bivalent vertex of colour 2 with an edge.
        for v in range(self.n_vertices, self.n_vertices + self.n_hairs + self.n_edges):
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
        return G

    def _hairy_to_bi_colored_hairy(self, G):
        # Translate hairy graphs to bi colored hairy graphs.
        # Generate all shuffles of the two hair types and assign the hairs to colour a and b respectively.
        hair_shuffles = ShuffleProduct(list(range(self.n_vertices, self.n_vertices + self.n_hairs_a)),
                                       list(range(self.n_vertices + self.n_hairs_a, self.n_vertices + self.n_hairs)))
        bi_colored_hairy_graphs = []
        vertices = list(range(0, self.n_vertices))
        for hairs in hair_shuffles:
            GG = copy(G)
            GG.relabel(vertices + hairs)
            bi_colored_hairy_graphs.append(GG)
        return bi_colored_hairy_graphs


class BiColoredHairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of bi colored hairy graph vector spaces with specified edge and hair parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_a_range (range): Range for the number of hairs of the first colour.
        - h_b_range (range): Range for the number of hairs of the second colour..
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs_a (bool): True for even hairs_a, False for odd hairs_a.
        - even_hairs_b (bool): True for even hairs_b, False for odd hairs_b.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_a_range, h_b_range, even_edges, even_hairs_a, even_hairs_b):
        """Initialize the sum vector space.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param h_a_range: Range for the number of hairs_a.
        :type h_a_range: range
        :param h_b_range: Range for the number of hairs_b.
        :type h_b_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: True for even hairs_a, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: True for even hairs_b, False for odd hairs_b.
        :type even_hairs_b: bool
        """
        self.v_range = v_range
        self.l_range = l_range
        self.h_a_range = h_a_range
        self.h_b_range = h_b_range
        self.sub_type = get_sub_type(even_edges, even_hairs_a, even_hairs_b)

        vs_list = []

        symmetry = not isinstance(
            self, GraphVectorSpace.DegSlice) and even_hairs_a == even_hairs_b

        for (v, l, h_a, h_b) in itertools.product(self.v_range, self.l_range, self.h_a_range, self.h_b_range):
            if symmetry and h_a < h_b:
                continue  # Symmetry between a and b hairs.
            vs_list.append(BiColoredHairyGraphVS(
                v, l, h_a, h_b, even_edges, even_hairs_a, even_hairs_b))
        super(BiColoredHairyGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('hairs_a', self.h_a_range),
                                   ('hairs_b', self.h_b_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s_%s" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


# ------- Operators --------
class ContractEdgesGO(HairyGraphComplex.ContractEdgesGO):
    """Contract edges graph operator.

    Operates on a hairy graph by contracting an edge not connected to a hair vertex and unifying the two adjacent vertices.

    Attributes:
        - sub_type (str): Graphs sub type of the domain.
    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: Domain vector space of the operator.
        :type domain: BiColoredHairyGraphVS
        :param target: Target vector space of the operator.
        :type target: BiColoredHairyGraphVS
        """
        self.sub_type = domain.sub_type
        super(ContractEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding contract edges graph operator.

        The contract edges operator reduces the number of vertices by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: BiColoredHairyGraphVS
        :param target: Potential target vector space of the operator.
        :type target: BiColoredHairyGraphVS
        :return: True if domain and target match to generate a corresponding contract edges graph operator.
        :rtype: bool
        """
        return domain.n_vertices - 1 == target.n_vertices and domain.n_loops == target.n_loops and \
            domain.n_hairs_a == target.n_hairs_a and domain.n_hairs_b == target.n_hairs_b \
            and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b):
        """Return a contract edges graph operator.

        :param n_vertices: Number of vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param n_hairs_a: Number of hairs_a.
        :type n_hairs_a: int
        :param n_hairs_b: Number of hairs_b.
        :type n_hairs_b: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: True for even hairs, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: True for even hairs, False for odd hairs_b.
        :type even_hairs_b: bool
        :return: Contract edges graph operator based on the specified domain vector space.
        :rtype: ContractEdgesGO
        """
        domain = BiColoredHairyGraphVS(
            n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b)
        target = BiColoredHairyGraphVS(n_vertices - 1, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a,
                                       even_hairs_b)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: BiColoredHairyGraphVS
        """
        super(ContractEdgesD, self).__init__(sum_vector_space,
                                             ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class SplitEdgesGO(GraphOperator.GraphOperator):
    """Split edges graph operator.

    Operate on a bi colored hairy graph by deleting an edge and adding a hair of colour a to one of the adjacent vertices and
    a hair of colour b to the other adjacent vertex.
    Only for graphs with odd edges, even hairs_a, and even hairs_b.

    Attributes:
        - sub_type (str): Graphs sub type of the domain.
    """

    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: BiColoredHairyGraphVS: Domain vector space of the operator.
        :param target: BiColoredHairyGraphVS: Target vector space of the operator.
        """
        self.sub_type = domain.sub_type
        super(SplitEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding contract edges graph operator.

        The split edges operator reduces the number of loops by one and increases the number of hairs by one for each
        colour of hairs.

        :param domain: BiColoredHairyGraphVS: Potential domain vector space of the operator.
        :param target: BiColoredHairyGraphVS: Potential target vector space of the operator.
        :return: bool: True if domain and target match to generate a corresponding split edges graph operator.
        """
        return domain.n_vertices == target.n_vertices and domain.n_loops - 1 == target.n_loops and \
            domain.n_hairs_a + 1 == target.n_hairs_a and domain.n_hairs_b + 1 == target.n_hairs_b \
            and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b):
        """Return a split edges graph operator.

        :param n_vertices: Number of vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param n_hairs_a: Number of hairs_a.
        :type n_hairs_a: int
        :param n_hairs_b: Number of hairs_b.
        :type n_hairs_b: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: True for even hairs, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: True for even hairs, False for odd hairs_b.
        :type even_hairs_b: bool
        :return: Split edges graph operator based on the specified domain vector space.
        :rtype: SplitEdgesGO
        """
        domain = BiColoredHairyGraphVS(
            n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b)
        target = BiColoredHairyGraphVS(n_vertices, n_loops - 1, n_hairs_a + 1, n_hairs_b + 1, even_edges, even_hairs_a,
                                       even_hairs_b)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "splitD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "splitD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
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
        return 'split edges'

    def operate_on(self, G):
        # Operate on a bi colored hairy graph by deleting an edge and adding a hair of colour a to one of the adjacent
        # vertices and a hair of colour b to the other adjacent vertex.
        sgn0 = -1 if G.order() % 2 else 1
        image = []
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # Only edges not connected to a hair-vertex can be split.
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices:
                continue
            G1 = copy(G)
            if not self.domain.even_edges:
                Shared.enumerate_edges(G1)
                e_label = G1.edge_label(u, v)
            G1.delete_edge((u, v))

            new_hair_idx_1 = self.domain.n_vertices + self.domain.n_hairs
            new_hair_idx_2 = new_hair_idx_1 + 1

            G1.add_vertex(new_hair_idx_1)
            G1.add_edge((u, new_hair_idx_1))
            G1.add_vertex(new_hair_idx_2)
            G1.add_edge((v, new_hair_idx_2))
            G2 = copy(G1)

            vertices = list(range(0, self.domain.n_vertices))
            vertices = [] if vertices is None else vertices
            start_idx_a = self.domain.n_vertices
            start_idx_b = self.domain.n_vertices + self.domain.n_hairs_a + 1
            hairs_a = list(range(start_idx_a + 1, start_idx_b))
            hairs_a = [] if hairs_a is None else hairs_a
            hairs_b = list(range(start_idx_b, new_hair_idx_2))
            hairs_b = [] if hairs_b is None else hairs_b
            p = vertices + hairs_a + hairs_b

            p1 = p + [start_idx_a, new_hair_idx_2]
            G1.relabel(p1)
            p2 = p + [new_hair_idx_2, start_idx_a]
            G2.relabel(p2)

            if not self.domain.even_edges:
                G1.set_edge_label(u, start_idx_a, e_label)
                G1.set_edge_label(v, new_hair_idx_2, G1.size() - 1)
                sgn1 = Shared.edge_perm_sign(G1)
                G2.set_edge_label(v, start_idx_a, e_label)
                G2.set_edge_label(u, new_hair_idx_2, G2.size() - 1)
                sgn2 = Shared.edge_perm_sign(G2)
            else:
                sgn1 = 1
                sgn2 = -1
            image.append((G1, sgn1*sgn0))
            image.append((G2, sgn2*sgn0))
        return image


class SplitEdgesD(GraphOperator.Differential):
    """Split edges differential.

    Only for graphs with odd edges, even hairs_a, and even hairs_b.
    """

    def __init__(self, sum_vector_space):
        """Initialize the split edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: BiColoredHairyGraphSumVS
        """
        super(SplitEdgesD, self).__init__(sum_vector_space,
                                          SplitEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'split edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_split_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_split_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


# ------- Graph Complex --------
class BiColoredHairyGC(GraphComplex.GraphComplex):
    """Graph complex for bi colored hairy graphs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_a_range (range): Range for the number of hairs of the first colour.
        - h_b_range (range): Range for the number of hairs of the second colour..
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs_a (bool): True for even hairs_a, False for odd hairs_a.
        - even_hairs_b (bool): True for even hairs_b, False for odd hairs_b.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_a_range, h_b_range, even_edges, even_hairs_a, even_hairs_b, differentials):
        """Initialize the graph complex.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param h_a_range: Range for the number of hairs_a.
        :type h_a_range: range
        :param h_b_range: Range for the number of hairs_b.
        :type h_b_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: True for even hairs_a, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: True for even hairs_b, False for odd hairs_b.
        :type even_hairs_b: bool
        :param differentials: List of differentials. Options: 'contract', 'split'.
        :type differentials: list(str)
        """
        self.v_range = v_range
        self.l_range = l_range
        self.h_a_range = h_a_range
        self.h_b_range = h_b_range
        self.sub_type = get_sub_type(even_edges, even_hairs_a, even_hairs_b)

        sum_vector_space = BiColoredHairyGraphSumVS(self.v_range, self.l_range, self.h_a_range, self.h_b_range,
                                                    even_edges, even_hairs_a, even_hairs_b)
        differential_list = []
        if not set(differentials).issubset(['contract', 'split']):
            raise ValueError(
                "Differentials for bi colored hairy graph complex: 'contract', 'split'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)
        if 'split' in differentials:
            split_edges_dif = SplitEdgesD(sum_vector_space)
            differential_list.append(split_edges_dif)
        super(BiColoredHairyGC, self).__init__(
            sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))
