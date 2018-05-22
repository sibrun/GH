
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


zero_hairs = False      # Option to include zero hairs in the hairy graph complexes.


# ------- Graph Vector Space --------
class BiColoredHairyGraphVS(GraphVectorSpace.GraphVectorSpace):

    def __init__(self, n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b):
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs_a = n_hairs_a
        self.n_hairs_b = n_hairs_b
        self.n_hairs = self.n_hairs_a + self.n_hairs_b
        self.even_edges = even_edges
        self.even_hairs_a = even_hairs_a
        self.even_hairs_b = even_hairs_b
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = get_sub_type(self.even_edges, self.even_hairs_a, self.even_hairs_b)
        super(BiColoredHairyGraphVS, self).__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(self.n_vertices + self.n_hairs, self.n_loops, self.even_edges)


    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and \
               self.n_hairs_a == other.n_hairs_a and self.n_hairs_b == other.n_hairs_b

    def get_basis_file_path(self):
        s = "gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops),
                                   ('hairs_a', self.n_hairs_a), ('hairs_b', self.n_hairs_b)])

    def get_partition(self):
        # all internal vertices are in color 1, the hair vertices are in color 2
        return [list(range(0, self.n_vertices)), list(range(self.n_vertices, self.n_vertices + self.n_hairs_a)),
                list(range(self.n_vertices + self.n_hairs_a, self.n_vertices + self.n_hairs))]

    def is_valid(self):
        # at least trivalent
        l = (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # all numbers positive
        l = l and self.n_vertices > 0 and self.n_loops >= 0 and \
            ((self.n_hairs_a >= 0 and self.n_hairs_b >= 0) if zero_hairs
             else (self.n_hairs_a > 0 and self.n_hairs_b > 0))
        # Can have at most a full graph
        l = l and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2
        # can have at most one hair per vertex
        l = l and self.n_vertices >= max(self.n_hairs_a, self.n_hairs_b)
        return l

    def get_work_estimate(self):
        if not self.is_valid():
            return 0
        # give estimate of number of graphs
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_generating_graphs(self):
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2
        # degree 1 pieces are hairs, degree 2 vertices are edges and are removed later
        # z switch prevents multiple hairs and multiple edges
        if not self.is_valid():
            return []
        n_vertices_1 = self.n_vertices
        n_vertices_2 = self.n_hairs + self.n_edges
        n_edges_bip = self.n_hairs + 2 * self.n_edges
        deg_range_1 = (3, n_edges_bip)
        deg_range_2 = (1, 2)
        bipartite_graphs = NautyInterface.list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2,
                                                                n_edges_bip)
        hairy_graphs = [self._bip_to_ordinary(G) for G in bipartite_graphs]
        list_of_lists = [self._hairy_to_bi_colored_hairy(G) for G in hairy_graphs]
        return list(itertools.chain.from_iterable(list_of_lists))

    def perm_sign(self, G, p):
        # the sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the hair-vertices
        sgn = self.ogvs.perm_sign(G, p)
        # compute the extra contribution from hairs if necessary
        if self.even_hairs_a == self.even_edges:
            hairs_a = p[self.n_vertices: self.n_vertices + self.even_hairs_a]
            if len(hairs_a) != 0:
                sgn *= Shared.Perm.shifted(hairs_a).signature()
        if self.even_hairs_b == self.even_edges:
            hairs_b = p[self.n_vertices + self.even_hairs_a : self.n_vertices + self.n_hairs]
            if len(hairs_b) != 0:
                sgn *= Shared.Perm.shifted(hairs_b).signature()
        return sgn

    def _bip_to_ordinary(self, G):
        # translates bipartite into ordinary graph by replacing a bivalent vertex of coulour 2 with an edge
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

    def _hairy_to_bi_colored_hairy(self, G):
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
    def __init__(self, v_range, l_range, h_a_range, h_b_range, even_edges, even_hairs_a, even_hairs_b):
        self.v_range = v_range
        self.l_range = l_range
        self.h_a_range = h_a_range
        self.h_b_range = h_b_range
        self.sub_type = get_sub_type(even_edges, even_hairs_a, even_hairs_b)

        vs_list = []
        for (v, l, h_a, h_b) in itertools.product(self.v_range, self.l_range, self.h_a_range, self.h_b_range):
            if even_hairs_a == even_hairs_b and h_a < h_b:
                continue # Symmetry between a and b hairs.
            vs_list.append(BiColoredHairyGraphVS(v, l, h_a, h_b, even_edges, even_hairs_a, even_hairs_b))
        super(BiColoredHairyGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('hairs_a', self.h_a_range),
                                   ('hairs_b', self.h_b_range)])


# ------- Operators --------
class ContractEdgesGO(HairyGraphComplex.ContractEdgesGO):
    def __init__(self, domain, target):
        self.sub_type = domain.sub_type
        super(ContractEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return domain.n_vertices - 1 == target.n_vertices and domain.n_loops == target.n_loops and \
               domain.n_hairs_a == target.n_hairs_a and domain.n_hairs_b == target.n_hairs_b \
               and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops,n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b):
        domain = BiColoredHairyGraphVS(n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b)
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
    def __init__(self, sum_vector_space):
        super(ContractEdgesD, self).__init__(sum_vector_space, ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.get_vs_list()[0].sub_type
        s = "cohomology_dim_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class SplitEdgesGO(GraphOperator.GraphOperator):
    def __init__(self, domain, target):
        self.sub_type = domain.sub_type
        super(SplitEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return domain.n_vertices == target.n_vertices and domain.n_loops - 1 == target.n_loops and \
               domain.n_hairs_a + 1 == target.n_hairs_a and domain.n_hairs_b + 1 == target.n_hairs_b \
               and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops,n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b):
        domain = BiColoredHairyGraphVS(n_vertices, n_loops, n_hairs_a, n_hairs_b, even_edges, even_hairs_a, even_hairs_b)
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

    def operate_on(self,G):
        sgn0 = -1 if G.order() % 2 else 1
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # only edges not connected to a hair-vertex can be cut
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
                G1.set_edge_label(v, new_hair_idx_2, G1.size() - 1 )
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
    def __init__(self, sum_vector_space):
        super(SplitEdgesD, self).__init__(sum_vector_space, SplitEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'split edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.get_vs_list()[0].sub_type
        s = "cohomology_dim_split_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


# ------- Graph Complex --------
class BiColoredHairyGC(GraphComplex.GraphComplex):
    def __init__(self, v_range, l_range, h_a_range, h_b_range, even_edges, even_hairs_a, even_hairs_b, differentials):
        self.v_range = v_range
        self.l_range = l_range
        self.h_a_range = h_a_range
        self.h_b_range = h_b_range
        self.sub_type = get_sub_type(even_edges, even_hairs_a, even_hairs_b)

        sum_vector_space = BiColoredHairyGraphSumVS(self.v_range, self.l_range, self.h_a_range, self.h_b_range,
                                                    even_edges, even_hairs_a, even_hairs_b)
        differential_list = []
        if not differentials <= {'contract', 'split'}:
            raise ValueError("Differentials for bi colored hairy graph complex: 'contract', 'split'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)
        if 'split' in differentials:
            split_edges_dif = SplitEdgesD(sum_vector_space)
            differential_list.append(split_edges_dif)
        super(BiColoredHairyGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))