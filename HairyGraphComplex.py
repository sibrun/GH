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


# ------- Hairy Graph Complex --------
class HairyGVS(GVS.GraphVectorSpace):

    def __init__(self, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs = n_hairs
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))
        super(HairyGVS,self).__init__()
        self.ogvs = OGC.OrdinaryGVS(self.n_vertices + self.n_hairs, self.n_loops, self.even_edges)

    def get_params(self):
        return (self.n_vertices, self.n_loops, self.n_hairs)

    def set_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % (self.n_vertices, self.n_loops, self.n_hairs)
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def set_plot_path(self):
        s = "gra%d_%d_%d" % (self.n_vertices, self.n_loops, self.n_hairs)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % (self.n_vertices, self.n_loops, self.n_hairs)
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def __str__(self):
        return "<Hairy graphs: %d vertices, %d loops, %d hairs, %s>" \
               % (self.n_vertices, self.n_loops, self.n_hairs, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs \
               and self.sub_type == other.sub_type

    def set_validity(self):
        # at least trivalent
        l = (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # all numbers positive
        l = l and self.n_vertices > 0 and self.n_loops >= 0 and self.n_hairs >= 0
        # Can have at most a full graph
        l = l and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2
        # can have at most one hair per vertex
        l = l and self.n_vertices >= self.n_hairs
        return l

    def set_partition(self):
        # all internal vertices are in color 1, the hair vertices are in color 2
        return [list(range(0, self.n_vertices)), list(range(self.n_vertices, self.n_vertices + self.n_hairs))]

    '''Produces a set of graphs whose isomorphism classes span the vector space. (Not necessarily freely!)'''
    def _generating_graphs(self):
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2
        # degree 1 pieces are hairs, degree 2 vertices are edges and are removed later
        # z switch prevents multiple hairs and multiple edges
        if not self.valid:
            return []
        n_vertices_1 = self.n_vertices
        n_vertices_2 = self.n_hairs + self.n_edges
        max_deg_1 = self.n_edges + 1
        n_edges_bip = self.n_hairs + 2 * self.n_edges
        L = NI.list_bipartite_graphs(n_vertices_1, n_vertices_2, (3, max_deg_1), (1,2), n_edges_bip)
        return [self._bip_to_ordinary(G) for G in L]

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

    def get_work_estimate(self):
        if not self.valid:
            return 0
        # give estimate of number of graphs
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)


# ------- Ordinary Graph Complex --------
class ContractDHairy(GO.GraphOperator):
    def __init__(self, domain, target):
        if domain.n_vertices != target.n_vertices + 1 or domain.n_loops != target.n_loops  or \
                        domain.n_hairs != target.n_hairs or domain.sub_type != target.sub_type:
            raise ValueError("Domain and target not consistent for contract edge operator")
        self.sub_type = sub_types.get((domain.even_edges, domain.even_hairs))
        super(ContractDHairy, self).__init__(domain, target)

    @classmethod
    def generate_operators(cls, vs_list):
        op_list = []
        for (domain, target) in itertools.product(vs_list, vs_list):
            if domain.n_vertices == target.n_vertices + 1 and domain.n_loops == target.n_loops \
                    and domain.n_hairs == target.n_hairs:
                op_list.append(cls(domain, target))
        return op_list

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        domain = HairyGVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs)
        target = HairyGVS(n_vertices - 1, n_loops, n_hairs, even_edges, even_hairs)
        return cls(domain, target)

    def get_params(self):
        return self.domain.get_params()

    def get_type(self):
        return 'contract edges'

    def set_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def set_rank_file_path(self):
        s = "contractD%d_%d_%d_rank.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d.txt.rank.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        if not self.valid:
            return 0
        return self.domain.n_edges * sqrt(self.target.get_dimension())

    def __str__(self):
        return "<Contract edges: domain: %s>" % str(self.domain)

    '''For G a graph returns a list of pairs (GG, x),
       such that (operator)(G) = sum x GG.'''
    def _operate_on(self,G):
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
        (v_range, l_range, h_range) = param_range
        return (range(min(v_range) + 1, max(v_range)), l_range, h_range)


# ------- Ordinary Graph Complex --------
class HairyGC(GC.GraphComplex):
    def __init__(self, v_range, l_range, h_range, even_edges, even_hairs):
        self.v_range = v_range
        self.l_range = l_range
        self.h_range = h_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = sub_types.get((self.even_edges, self.even_hairs))

        vs_list = [HairyGVS(v, l, h, self.even_edges, self.even_hairs) for
                   (v, l, h) in itertools.product(self.v_range, self.l_range, self.h_range)]
        op_list = ContractDHairy.generate_operators(vs_list)
        super(HairyGC, self).__init__(vs_list, op_list)

    def get_type(self):
        e = 'even edges' if self.even_edges else 'odd edges'
        h = 'even hairs' if self.even_hairs else 'odd hairs'
        return e + ', ' + h

    def get_params_range(self):
        return [self.v_range, self.l_range, self.h_range]

    def get_params_names(self):
        return ('vertices', 'loops', 'hairs')

    def __str__(self):
        return "<Hairy graph complex with %s and parameter range: vertices: %s, loops: %s, hairs: %s>" \
               % (self.sub_type, str(self.v_range), str(self.l_range), str(self.h_range))

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
            dim_dict.update({(vs.n_vertices, vs.n_loops, vs.n_hairs): cohomology_dim.get(vs)})
        param_range = ContractDHairy.transform_param_range((self.v_range, self.l_range, self.h_range))
        return (dim_dict, param_range)

    def plot_cohomology_dim(self):
        (dim_dict, param_range) = self.get_cohomology_dim()
        (v_range, l_range, h_range) = param_range
        path = self.get_cohomology_plot_path()
        Display.plot_3d_array(dim_dict, 'vertices', v_range, 'loops', l_range, 'hairs', h_range, path)
