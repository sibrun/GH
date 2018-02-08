import itertools
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH
import NautyInterface as NI
import Display
import Parameters


graph_type = "ordinary"
sub_types = {True: "even_edges", False: "odd_edges"}


# ------- Ordinary Graph Vector Space --------
class OrdinaryGVS(GVS.GraphVectorSpace):

    def __init__(self, n_vertices, n_loops, even_edges):
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        super(OrdinaryGVS,self).__init__()

    def get_params(self):
        return (self.n_vertices, self.n_loops)

    def set_basis_file_path(self):
        s = "gra%d_%d.g6" % (self.n_vertices, self.n_loops)
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def set_plot_path(self):
        s = "gra%d_%d" % (self.n_vertices, self.n_loops)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d.g6" % (self.n_vertices, self.n_loops)
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def set_validity(self):
        return (3 * self.n_vertices <= 2 * self.n_edges) and self.n_vertices > 0 and self.n_loops >= 0 \
               and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2

    def set_partition(self):
        return None

    def get_work_estimate(self):
        if not self.valid:
            return 0
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def __str__(self):
        return "<Ordinary graphs: %d vertices, %d loops, %s>" % (self.n_vertices, self.n_loops, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops \
               and self.even_edges == other.even_edges

    def _generating_graphs(self):
        if not self.valid:
            return []
        return NI.list_simple_graphs(self.n_vertices, self.n_edges)

    def perm_sign(self, G, p):
        if self.even_edges:
            # The sign is (induced sign on vertices) * (induced sign edge orientations)
            sign = SH.Perm(p).sign()
            for (u, v) in G.edges(labels=False):
                # we assume the edge is always directed from the larger to smaller index
                if (u < v and p[u] > p[v]) or (u > v and p[u] < p[v]):
                    sign *= -1
            return sign
        else:
            # The sign is (induced sign of the edge permutation)
            # we assume the edges are always lex ordered
            # for the computation we use that G.edges() returns the edges in lex ordering
            # we first label the edges on a copy of G lexicographically
            G1 = copy(G)
            for (j,e) in enumerate(G1.edges(labels=False)):
                (u, v) = e
                G1.set_edge_label(u, v, j)

            # we permute the graph, and read of the new labels
            G1.relabel(p, inplace=True)
            return SH.Perm([j for (u, v, j) in G1.edges()]).sign()


# ------- Contraction Operator --------
class ContractDOrdinary(GO.GraphOperator):

    def __init__(self, domain, target):
        if domain.n_vertices != target.n_vertices+1 or domain.n_loops != target.n_loops \
                or domain.even_edges != target.even_edges:
            raise ValueError("Domain and target not consistent for contract edge operator")
        self.sub_type = sub_types.get(domain.even_edges)
        super(ContractDOrdinary, self).__init__(domain, target)

    @classmethod
    def generate_operators(cls, vs_list):
        op_list = []
        for (domain, target) in itertools.product(vs_list, vs_list):
            if domain.n_vertices == target.n_vertices + 1 and domain.n_loops == target.n_loops:
                op_list.append(cls(domain, target))
        return op_list

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, even_edges):
        domain = OrdinaryGVS(n_vertices, n_loops, even_edges)
        target = OrdinaryGVS(n_vertices - 1, n_loops, even_edges)
        return cls(domain, target)

    def get_params(self):
        return self.domain.get_params()

    def get_type(self):
        return 'contract edges'

    def set_matrix_file_path(self):
        s = "contractD%d_%d.txt" % (self.domain.n_vertices, self.domain.n_loops)
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def set_rank_file_path(self):
        s = "contractD%d_%d_rank.txt" % (self.domain.n_vertices, self.domain.n_loops)
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d.txt" % (self.domain.n_vertices, self.domain.n_loops)
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d.txt.rank.txt" % (self.domain.n_vertices, self.domain.n_loops)
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        if not self.valid:
            return 0
        return self.domain.n_edges * sqrt(self.target.get_dimension())

    def __str__(self):
        return "<Contract edges: domain: %s>" % str(self.domain)

    def _operate_on(self,G):
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            r = range(0,self.domain.n_vertices)
            p = list(r)
            p[0] = u
            p[1] = v
            idx = 2
            for j in r:
                if j == u or j== v:
                    continue
                else:
                    p[idx] = j
                    idx +=1

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

    @staticmethod
    def transform_param_range(param_range):
        (v_range, l_range) = param_range
        return (range(min(v_range) + 1, max(v_range)), l_range)


# ------- Ordinary Graph Complex --------
class OrdinaryGC(GC.GraphComplex):
    def __init__(self, v_range, l_range, even_edges):
        self.v_range = v_range
        self.l_range = l_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        vs_list = [OrdinaryGVS(v, l, self.even_edges) for (v, l) in itertools.product(self.v_range, self.l_range)]
        op_list = ContractDOrdinary.generate_operators(vs_list)
        super(OrdinaryGC, self).__init__(vs_list, op_list)

    def get_type(self):
        return 'even edges' if self.even_edges else 'odd edges'

    def get_params_range(self):
        return (self.v_range, self.l_range)

    def get_params_names(self):
        return('vertices', 'loops')

    def __str__(self):
        return "<Ordinary graph complex with %s and parameter range: vertices: %s, loops: %s>" \
               % (self.sub_type, str(self.v_range), str(self.l_range))

    def set_info_file_path(self):
        s = "graph_complex.txt"
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_cohomology_plot_path(self):
        s = "cohomology_dim_%s_%s.png" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

    def get_cohomology_dim(self):
        cohomology_dim = self.get_general_cohomology_dim_dict()
        dim_dict = dict()
        for vs in self.vs_list:
            dim_dict.update({(vs.n_vertices, vs.n_loops): cohomology_dim.get(vs)})
        param_range = ContractDOrdinary.transform_param_range((self.v_range, self.l_range))
        return(dim_dict, param_range)

    def plot_cohomology_dim(self):
        (dim_dict, param_range) = self.get_cohomology_dim()
        (v_range, l_range) = param_range
        path = self.get_cohomology_plot_path()
        Display.plot_2d_array(dim_dict, 'vertices', v_range, 'loops', l_range, path)
