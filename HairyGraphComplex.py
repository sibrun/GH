import itertools
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH
import NautyInterface as NI
import StoreLoad as SL
import Display
import OrdinaryGraphComplex as OGC


data_dir = "data"
plots_dir = "plots"
ref_data_dir = "data_ref"
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

    def set_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % (self.n_vertices, self.n_loops, self.n_hairs)
        return os.path.join(data_dir, graph_type, self.sub_type, s)

    def set_plot_path(self):
        s = "cohomology%d_%d_%d.png" % (self.n_vertices, self.n_loops, self.n_hairs)
        return os.path.join(plots_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d_%d.g6" % (self.n_vertices, self.n_loops, self.n_hairs)
        return os.path.join(ref_data_dir, graph_type, self.sub_type, s)

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
        # give estimate of number of graphs
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)


# ------- Ordinary Graph Complex --------
class ContractGO(GO.GraphOperator):
    def __init__(self, domain, target):
        if domain.n_vertices != target.n_vertices+1 or domain.n_loops != target.n_loops \
                or domain.sub_type != target.sub_type:
            raise ValueError("Domain and target not consistent for contract edge operator")
        self.sub_type = sub_types.get((domain.even_edges, domain.even_hairs))
        super(ContractGO, self).__init__(domain, target)

    @classmethod
    def get_operators(cls, vs_list):
        op_list = []
        for (domain, target) in itertools.product(vs_list, vs_list):
            if domain.n_vertices == target.n_vertices + 1 and domain.n_loops == target.n_loops:
                op_list.append(cls(domain, target))
        return op_list

    @classmethod
    def get_operator(cls, n_vertices, n_loops, n_hairs, even_edges, even_hairs):
        domain = HairyGVS(n_vertices, n_loops, n_hairs, even_edges, even_hairs)
        target = HairyGVS(n_vertices - 1, n_loops, n_hairs, even_edges, even_hairs)
        return cls(domain, target)

    def set_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(data_dir, graph_type, self.sub_type, s)

    def set_rank_file_path(self):
        s = "contractD%d_%d_%d_rank.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d.txt.rank.txt" % (self.domain.n_vertices, self.domain.n_loops, self.domain.n_hairs)
        return os.path.join(ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
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
        op_list = ContractGO.get_operators(vs_list)
        super(HairyGC, self).__init__(vs_list, op_list)

    def __str__(self):
        return "<Hairy graph complex with %s and parameter range: vertices: %s, loops: %s>" \
               % (self.sub_type, str(self.v_range), str(self.l_range))

    def _set_info_file_path(self):
        s = "graph_complex.txt"
        return os.path.join(data_dir, graph_type, self.sub_type, s)

    def get_cohomology_plot_path(self):
        s = "cohomology_dim_%s_%s.png" % (graph_type, self.sub_type)
        return os.path.join(plots_dir, graph_type, self.sub_type, s)

    def get_cohomology_file_path(self):
        s = "cohomology_dim_%s_%s.p" % (graph_type, self.sub_type)
        return os.path.join(data_dir, graph_type, self.sub_type, s)

    def compute_cohomology_dim(self):
        self._compute_cohomology_dim()
        dim_dict = dict()
        for vs in self.vs_list:
            dim_dict.update({(vs.n_vertices, vs.n_loops, vs.n_hairs): self.cohomology_dim.get(vs)})
        path = self.get_cohomology_file_path()
        v_range = range(min(self.v_range)+1,max(self.v_range))
        SL.pickle_store((dim_dict, v_range, self.l_range, self.h_range), path)

    def get_cohomology_dim(self):
        if not self.exists_cohomology_file():
            raise SL.NotBuiltError("Cannot load cohomology dimensions, No cohomology file found for %s: " % str(self))
        (dim_dict, v_range, l_range, h_range) = SL.pickle_load(self.get_cohomology_file_path())
        return dim_dict

    def plot_cohomology_dim(self):
        if not self.exists_cohomology_file():
            raise SL.NotBuiltError("Cannot load cohomology dimensions, No cohomology file found for %s: " % str(self))
        (dim_dict, v_range, l_range, h_range) = SL.pickle_load(self.get_cohomology_file_path())
        path = self.get_cohomology_plot_path()
        Display.plot_3d_array(dim_dict, 'vertices', v_range, 'loops', l_range, 'hairs', h_range, path)
