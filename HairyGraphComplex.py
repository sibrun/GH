from abc import ABCMeta, abstractmethod
import math
import itertools
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH
import NautyInterface as NI
import Display
import OrdinaryGraphComplex as OGC
import Parameters


graph_type = "hairy"
sub_types = {(True, True): "even_edges_even_hairs", (True, False): "even_edges_odd_hairs",
             (False, True): "odd_edges_even_hairs", (False, False): "odd_edges_odd_hairs"}


# ------- Graph Vector Space --------
class HairySubGVS(GVS.SubGraphVectorSpace):

    def __init__(self, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs = n_hairs
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))
        super(HairySubGVS,self).__init__()
        self.ogvs = OGC.OrdinarySubGVS(self.n_vertices + self.n_hairs, self.n_loops, self.even_edges)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs \
               and self.sub_type == other.sub_type

    def get_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % self.get_params_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_plot_path(self):
        s = "gra%d_%d_%d" % self.get_params_tuple()
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % self.get_params_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_params_dict(self):
        return {'vertices': self.n_vertices, 'loops': self.n_loops, 'hairs': self.n_hairs}

    def get_params_tuple(self):
        return (self.n_vertices, self.n_loops, self.n_hairs)

    def get_partition(self):
        # all internal vertices are in color 1, the hair vertices are in color 2
        return [list(range(0, self.n_vertices)), list(range(self.n_vertices, self.n_vertices + self.n_hairs))]

    def is_valid(self):
        # at least trivalent
        l = (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # all numbers positive
        l = l and self.n_vertices > 0 and self.n_loops >= 0 and self.n_hairs >= 0
        # Can have at most a full graph
        l = l and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2
        # can have at most one hair per vertex
        l = l and self.n_vertices >= self.n_hairs
        return l

    def get_work_estimate(self):
        if not self.is_valid():
            return 0
        # give estimate of number of graphs
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    '''Produces a set of graphs whose isomorphism classes span the vector space. (Not necessarily freely!)'''
    def get_generating_graphs(self):
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2
        # degree 1 pieces are hairs, degree 2 vertices are edges and are removed later
        # z switch prevents multiple hairs and multiple edges
        if not self.is_valid():
            return []
        n_vertices_1 = self.n_vertices
        n_vertices_2 = self.n_hairs + self.n_edges
        max_deg_1 = self.n_edges + 1
        n_edges_bip = self.n_hairs + 2 * self.n_edges
        L = NI.list_bipartite_graphs(n_vertices_1, n_vertices_2, (3, max_deg_1), (1,2), n_edges_bip)
        return [self._bip_to_ordinary(G) for G in L]

    '''For G a graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
       Here vertex j becomes vertex p[j] in the new graph.'''
    def perm_sign(self, G, p):
        # the sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the hair-vertices
        sgn = self.ogvs.perm_sign(G, p)
        # compute the extra contribution from hairs if necessary
        if self.even_hairs == self.even_edges:
            hairs = p[self.n_vertices:]
            sgn *= SH.Perm.shifted(hairs).sign()
        return sgn

    def _bip_to_ordinary(self, G):
        #return G
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


class HairyGVS(GVS.GraphVectorSpace):
    def __init__(self, v_range, l_range, h_range, even_edges, even_hairs):
        self.v_range = v_range
        self.l_range = l_range
        self.h_range = h_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))

        vs_list = [HairySubGVS(v, l, h, self.even_edges, self.even_hairs) for
                   (v, l, h) in itertools.product(self.v_range, self.l_range, self.h_range)]
        super(HairyGVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_params_range_dict(self):
        return {'vertices': self.v_range, 'loops': self.l_range, 'hairs': self.h_range}


# ------- Gradings --------
class VertexGrading(GVS.Grading):
    def __init__(self, graph_vector_space):
        super(VertexGrading, self).__init__(graph_vector_space)

    def get_deg_idx(self, ordinary_sub_gvs):
        return ordinary_sub_gvs.n_vertices


class HairGrading(GVS.Grading):
    def __init__(self, graph_vector_space):
        super(HairGrading, self).__init__(graph_vector_space)

    def get_deg_idx(self, ordinary_sub_gvs):
        return ordinary_sub_gvs.n_hairs


# ------- Operators --------
class ContractEdgesGO(GO.GraphOperator):
    def __init__(self, domain, target):
        if not ContractEdgesGO.is_match(domain, target):
            raise ValueError("Domain and target not consistent for contract edge operator")
        self.sub_type = domain.sub_type
        super(ContractEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return domain.n_vertices == target.n_vertices + 1 and domain.n_loops == target.n_loops \
               and domain.n_hairs == target.n_hairs and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        domain = HairySubGVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs)
        target = HairySubGVS(n_vertices - 1, n_loops, n_hairs, even_edges, even_hairs)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % self.domain.get_params_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_rank.txt" % self.domain.get_params_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % self.domain.get_params_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d.txt.rank.txt" % self.domain.get_params_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        if not self.is_valid():
            return 0
        target_dim = self.target.get_dimension()
        if target_dim == 0:
            return 0
        return self.domain.n_edges * math.log(self.target.get_dimension(), 2)

    def get_type(self):
        return 'contract edges'

    '''For G a graph returns a list of pairs (GG, x),
       such that (operator)(G) = sum x GG.'''
    def operate_on(self,G):
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # only edges not connected to a hair-vertex can be contracted
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices:
                continue
            r = range(0,self.domain.n_vertices + self.domain.n_hairs)
            p = list(r)
            p[0] = u
            p[1] = v
            idx = 2
            for j in r:
                if j == u or j == v:
                    continue
                else:
                    p[idx] = j
                    idx += 1

            pp = SH.Perm(p).inverse()
            sgn = self.domain.perm_sign(G, pp)
            G1 = copy(G)
            G1.relabel(pp, inplace=True)

            for (j, ee) in enumerate(G1.edges(labels=False)):
                a, b = ee
                G1.set_edge_label(a,b,j)
            previous_size = G1.size()
            G1.merge_vertices([0,1])
            if (previous_size - G1.size()) != 1:
                continue
            G1.relabel(list(range(0,G1.order())), inplace=True)
            if not self.domain.even_edges:
                p = [j for (a, b, j) in G1.edges()]
                sgn *= Permutation(p).signature()
            image.append((G1, sgn))
        return image


class ContractEdgesD(GO.Differential):
    def __init__(self, vector_space):
        super(ContractEdgesD, self).__init__(vector_space, ContractEdgesGO.generate_op_matrix_list(vector_space))

    def get_type(self):
        return 'contract edges'


class MergeHairsGO(GO.GraphOperator):
    def __init__(self, domain, target):
        if not MergeHairsGO.is_match(domain, target):
            raise ValueError("Domain and target not consistent for contract edge operator")
        self.sub_type = domain.sub_type
        super(MergeHairsGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return domain.n_vertices == target.n_vertices and domain.n_loops == target.n_loops + 1 \
               and domain.n_hairs == target.n_hairs - 1 and domain.sub_type == target.sub_type

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        domain = HairySubGVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs)
        target = HairySubGVS(n_vertices, n_loops - 1, n_hairs + 1, even_edges, even_hairs)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "splitD%d_%d_%d.txt" % self.domain.get_params_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "splitD%d_%d_%d_rank.txt" % self.domain.get_params_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        pass

    def get_ref_rank_file_path(self):
        pass

    def get_work_estimate(self):
        if not self.is_valid():
            return 0
        target_dim = self.target.get_dimension()
        if target_dim == 0:
            return 0
        return self.domain.n_edges * math.log(self.target.get_dimension(), 2)

    def get_type(self):
        return 'split edges'

    '''For G a graph returns a list of pairs (GG, x),
       such that (operator)(G) = sum x GG.'''
    def operate_on(self,G):
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # only edges not connected to a hair-vertex can be split
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices:
                continue
            r = range(0,self.domain.n_vertices + self.domain.n_hairs)
            p = list(r)
            p[0] = u
            p[1] = v
            idx = 2
            for j in r:
                if j == u or j == v:
                    continue
                else:
                    p[idx] = j
                    idx += 1

            pp = SH.Perm(p).inverse()
            sgn = self.domain.perm_sign(G, pp)
            G1 = copy(G)
            G1.relabel(pp, inplace=True)

            for (j, ee) in enumerate(G1.edges(labels=False)):
                a, b = ee
                G1.set_edge_label(a,b,j)
            previous_size = G1.size()
            G1.merge_vertices([0,1])
            if (previous_size - G1.size()) != 1:
                continue
            G1.relabel(list(range(0,G1.order())), inplace=True)
            if not self.domain.even_edges:
                p = [j for (a, b, j) in G1.edges()]
                sgn *= Permutation(p).signature()
            image.append((G1, sgn))
        return image


class SplitEdgesD(GO.Differential):
    def __init__(self, vector_space):
        super(SplitEdgesD, self).__init__(vector_space, SplitEdgesGO.generate_op_matrix_list(vector_space))

    def get_type(self):
        return 'split edges'


# ------- Graph Complexes --------
class HairyContractEdgesGC(GC.GraphComplex):
    def __init__(self, v_range, l_range, h_range, even_edges, even_hairs):
        self.v_range = v_range
        self.l_range = l_range
        self.h_range = h_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))

        vector_space = HairyGVS(v_range, l_range, h_range, even_edges, even_hairs)
        grading = VertexGrading(vector_space)
        differential = ContractEdgesD(vector_space)
        super(HairyContractEdgesGC, self).__init__(vector_space, grading, differential)

    def get_cohomology_plot_path(self):
        s = "cohomology_dim_%s_%s.png" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

    def plot_cohomology_dim(self):
        dim_dict = self.differential.get_cohomology_dim()
        plot_path = self.get_cohomology_plot_path()
        param_labels_ranges = ('vertices', self.v_range, 'loops', self.l_range, 'hairs', self.h_range)
        Display.plot_3d_array(dim_dict, param_labels_ranges, plot_path)

