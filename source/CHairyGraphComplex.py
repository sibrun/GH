"""Graph complexes based on simple graphs with numbered (colored) hairs.
Each hair has its own color.
These graphs compute W_0 H_c(M_g,n), with n the number of hairs, g the loop order.
Implemented Differentials: Contract edges.
The first vertices correspond to internal vertices, and the last to the hairs.
"""


__all__ = ['graph_type', 'sub_types', 'CHairyGraphVS', 'CHairyGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
           'CHairyGC']

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
import SymmetricGraphComplex
import GCDimensions

graph_type = "chairy"

# Option to include zero hairs in the hairy graph complexes.
zero_hairs = False


def dict_to_list(d, n):
    return [(d[j] if j in d else j) for j in range(n)]


# ------- Graph Vector Space --------
class CHairyGraphVS(SymmetricGraphComplex.SymmetricGraphVectorSpace):
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
        - sub_type (str): Sub type of graphs.
        - ogvs (OrdinaryGraphComplex.OrdinaryGVS): Ordinary graph vector space without hairs.

    """

    def __init__(self, n_vertices, n_loops, n_hairs, even_edges):
        """Initialize the hairy graph vector space.

        :param n_vertices: Number of internal vertices. TODO: sure it is only internal?
        :type n_vertices: int
        :param n_loops: the genus of the graph.
        :type n_loops: int
        :param n_hairs: Number of hairs. They are distinguishable, numbered 1,...,n
        :type n_hairs: int
        :param n_ws: Number of w decorated hairs.
        :type even_edges: int
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs = n_hairs
        self.even_edges = even_edges
        self.sub_type = "even_edges" if even_edges else "odd_edges"

        # we count only the internal edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        super(CHairyGraphVS, self).__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(
            self.n_vertices + self.n_hairs, self.n_loops, even_edges)

    def get_type(self):
        return 'chgraphs'

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs and self.even_edges == other.even_edges

    def __hash__(self):
        return hash("chgra%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "chgra%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "chgra%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('hairs', self.n_hairs)])

    def get_partition(self):
        # All internal vertices are in color 1, the single eps vertex in color 2, the w vertex in color 3
        # and the hair vertices are in colors 4,...,n+3.
        return [list(range(0, self.n_vertices))] + [[j] for j in range(self.n_vertices, self.n_vertices + self.n_hairs)]

    def plot_graph(self, G):
        GG = Graph(G)  # , loops=True)

        return GG.plot(partition=self.get_partition(), vertex_labels=True)

    def is_valid(self):
        # At least trivalent internal vertices.
        l = (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # Nonnegative number of vertices, non negative number of loops, non-negative or positive number of hairs,
        l = l and self.n_vertices >= 0 and self.n_loops >= 0 and self.n_hairs >= 0
        # At most a full graph.
        l = l and self.n_edges <= (self.n_vertices) * (self.n_vertices - 1) / 2
        return l

    def get_work_estimate(self):
        # TODO
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return GCDimensions.get_chairy_dim_estimate(self.n_vertices, self.n_loops, self.n_hairs)
        # return (self.n_vertices ** self.n_hairs) * binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_hairy_graphs(self, nvertices, nloops, nhairs, include_novertgraph=false):
        """ Produces all connected hairy graphs with nhairs hairs, that are the last vertices in the ordering.
        Graphs can have multiple hairs, but not tadpoles or multiple edges.
        :param include_novertgraph: Whether to include the graph with one edge and no vertices as a two-hair graph
        :type include_novertgraph: bool
        """
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2.
        # Degree 1 vertices are hairs, degree 2 vertices are edges and are removed later.
        nedges = nloops + nvertices - 1  # number of internal edges
        n_vertices_1 = nvertices
        n_vertices_2 = nhairs + nedges
        n_edges_bip = nhairs + 2 * nedges
        deg_range_1 = (3, n_edges_bip + 1)
        deg_range_2 = (1, 2)

        # check if valid
        unordered = []
        if (nvertices >= 1 and nloops >= 0 and nhairs >= 0 and n_edges_bip >= n_vertices_2
            and n_edges_bip <= 2*n_vertices_2 and n_edges_bip >= 3 * n_vertices_1
                and n_edges_bip <= n_vertices_1 * n_vertices_2):
            bipartite_graphs = NautyInterface.list_bipartite_graphs2(
                n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip)
            unordered = [self._bip_to_ordinary(
                G, nvertices, nedges, nhairs) for G in bipartite_graphs]
        # Produce all permutations of the hairs
        # all_perm = [ range(0,nvertices) + p for p in Permutations(range(nvertices, nvertices+nhairs)) ]
        # return [G.relabel(p, inplace=False) for p in all_perm ]
        if include_novertgraph and nvertices == 0 and nhairs == 2 and nloops == 0:
            unordered.append(Graph([(0, 1)]))
        return unordered

    def get_generating_graphs(self):
        # The routines above produce all hairy graphs, we just have to permute the hair labels

        # Produce all permutations of the hairs
        all_perm = [list(range(0, self.n_vertices)) + list(p)
                    for p in itertools.permutations(range(self.n_vertices, self.n_vertices+self.n_hairs))]

        return (G.relabel(p, inplace=False) for G in self.get_hairy_graphs(self.n_vertices, self.n_loops, self.n_hairs)
                for p in all_perm)

    def perm_sign(self, G, p):
        # The sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the hair-vertices.
        sgn = self.ogvs.perm_sign(G, p)

        # Compute the extra contribution from hairs.
        # if self.even_hairs == self.even_edges:
        #     hairs = p[self.n_vertices:]
        #     if len(hairs) != 0:
        #         sgn *= Shared.Perm.shifted(hairs).signature()
        return sgn

    def _bip_to_ordinary(self, G, nvertices, nedges, nhairs):
        # Translates bipartite into ordinary graph by replacing a bivalent vertex of colour 2 with an edge.
        for v in range(nvertices, nvertices + nhairs + nedges):
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
        G.relabel(range(0, G.order()))
        return G

    def get_n(self):
        return self.n_hairs

    def vertex_permutation_from_permutation(self, p):
        return list(range(0, self.n_vertices)) + [j+self.n_vertices-1 for j in p]

    def get_isotypical_projector(self, rep_index):
        return SymmProjector(self, rep_index)


class CHairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of hairy graph vector spaces with specified number of omega hairs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_range (range): Range for the number of hairs.
        - w_range (range): number of omega hairs
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_range, even_edges):
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
        self.sub_type = "even_edges" if even_edges else "odd_edges"

        vs_list = [CHairyGraphVS(v, l, h, self.even_edges) for
                   (v, l, h) in itertools.product(self.v_range, self.l_range, self.h_range)]
        super(CHairyGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return 'chairy graphs'

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('hairs', self.h_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s" % graph_type
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


# ------- Operators --------
class ContractEdgesGO(SymmetricGraphComplex.SymmetricGraphOperator):
    """Contract edges graph operator.

    Operate on a w-hairy graph by contracting an edge not connected to a hair vertex and unifying the two adjacent vertices.

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
            and domain.n_hairs == target.n_hairs

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, even_edges):
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
        domain = CHairyGraphVS(n_vertices, n_loops, n_hairs, even_edges)
        target = CHairyGraphVS(n_vertices - 1, n_loops, n_hairs, even_edges)
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
            (domain_dim, dimtarget_dim) = (
                self.domain.get_dimension(), self.target.get_dimension())
        except StoreLoad.FileNotFoundError:
            return 0
        if domain_dim == 0 or dimtarget_dim == 0:
            return 0
        return self.domain.n_edges * domain_dim * math.log(self.target.get_dimension(), 2)

    def get_type(self):
        return 'contract edges'

    def operate_on(self, G):
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image = []
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # only edges not connected to a hair-vertex can be contracted
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices:
                continue
            pp = Shared.permute_to_left((u, v), range(
                0, self.domain.n_vertices + self.domain.n_hairs))
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

    def restrict_to_isotypical_component(self, rep_index):
        #opP = self.domain.get_isotypical_projector(rep_index)
        return RestrictedContractEdgesGO(self, rep_index)


class RestrictedContractEdgesGO(SymmetricGraphComplex.SymmetricRestrictedOperatorMatrix):
    # def __init__(opD, opP):

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_r%d.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_r%d_rank.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_work_estimate(self):
        return self.opD.get_work_estimate()

    def is_match(self, domain, target):
        return ContractEdgesGO.is_match(domain.vs, target.vs) and domain.rep_index == target.rep_index


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: HairyGraphSumVS
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
        s = "info_contract_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class RestrictedContractEdgesD(SymmetricGraphComplex.SymmetricDifferential):
    def get_type(self):
        return 'isotypical contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "info_contract_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class SymmProjector(SymmetricGraphComplex.SymmetricProjectionOperator):
    """This class encodes the projector to an isotypical component of the symmetric group action
        by permuting numbered hairs.
        Warning: The matrix stores not the projector, but projector * n_hairs! / rep_dimension??, to have integral matrices.

    Attributes:
        - sub_type(str): Graphs sub type of the domain.
    """

    # def norm_permutation(self, p):
    #     """Returns the permutation on the vertices of a graph corresponding to a permutation of letters 1,...,n.
    #     :param p: a permutation
    #     """
    #     nn = sum(p)
    #     return list(range(0, self.domain.n_vertices)) + [j+self.domain.n_vertices-1 for j in p]

    def __init__(self, domain, rep_index):
        """Initialize the domain and target vector space of the contract edges graph operator.

        : param domain: Domain vector space of the operator.
        : type domain: HairyGraphVS
        : param rep_index: The index of the representation in the list produced by Partitions(h).
        : type rep_index: int
        """
        self.sub_type = domain.sub_type

        super(SymmProjector, self).__init__(domain, rep_index)

        # fill in representation and character
        # nn = domain.n_hairs
        # self.rep_partition = Partitions(nn)[rep_index]
        # self.norm_char_perm = [(symmetrica.charvalue(self.rep_partition, p.cycle_type(
        # )), self.norm_permutation(p)) for p in Permutations(nn)]

        # print(self.norm_char_perm)
    @staticmethod 
    def generate_operator(n_vertices, n_loops, n_hairs, even_edges, rep_index):
        return CHairyGraphVS(n_vertices, n_loops, n_hairs, even_edges).get_isotypical_projector(rep_index)
    
    # @staticmethod
    # def is_match(domain, target):
    #     """Check whether domain and target match to generate a corresponding contract edges graph operator.

    #     The contract edges operator reduces the number of vertices by one.

    #     : param domain: Potential domain vector space of the operator.
    #     : type domain: HairyGraphVS
    #     : param target: Potential target vector space of the operator.
    #     : type target: HairyGraphVS
    #     : return: True if domain and target match to generate a corresponding contract edges graph operator.
    #     : rtype: bool
    #     """
    #     return domain == target

    def get_ordered_param_dict2(self):
        do = self.domain
        return Shared.OrderedDict([('vertices', do.n_vertices), ('loops', do.n_loops), ('hairs', do.n_hairs), ('rep_index', self.rep_index)])

    def get_matrix_file_path(self):
        s = "projectionO%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "projectionO%d_%d_%d_%d_rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "projectionO%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "projectionO%d_%d_%d_%d.txt.rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    # def get_type(self):
    #     return 'projection operator'

    # def operate_on(self, G):
    #     # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
    #     image = []
    #     for (c, p) in self.norm_char_perm:
    #         # c is char value, p is permutation
    #         G1 = copy(G)
    #         sgn = self.domain.ogvs.perm_sign(G1, p)
    #         G1.relabel(p, inplace=True)
    #         image.append((G1, sgn * c))

    #     return image


# ------- Graph Complex --------
class CHairyGC(GraphComplex.GraphComplex):
    """Graph complex for hairy graphs.

    Attributes:
        - v_range(range): Range for the number of vertices.
        - l_range(range): Range for the number of loops.
        - h_range(range): Range for the number of hairs.
        - w_range(range): Range for the number of omega hairs
        - sub_type(str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_range, even_edges, differentials):
        """Initialize the graph complex.

        : param v_range: Range for the number of vertices.
        : type v_range: range
        : param l_range: Range for the number of loops.
        : type l_range: range
        : param h_range: Range for the number of hairs.
        : type  h_range: range
        : param even_edges: True for even edges, False for odd edges.
        : type even_edges: bool
        : param even_hairs: True for even hairs, False for odd hairs.
        : type even_hairs: bool
        : param differentials: List of differentials. Options: 'contract', 'et1h'.
        : type differentials: list(str)
        """
        self.v_range = v_range
        self.l_range = l_range
        self.h_range = h_range
        self.even_edges = even_edges
        self.sub_type = "even_edges" if even_edges else "odd_edges"

        sum_vector_space = CHairyGraphSumVS(
            self.v_range, self.l_range, self.h_range, even_edges)
        differential_list = []
        if not set(differentials).issubset(['contract', 'contract_iso']):
            raise ValueError(
                "Differentials for hairy graph complex: 'contract', 'contract_iso'")
        contract_edges_dif = ContractEdgesD(sum_vector_space)
        if 'contract' in differentials:
            differential_list.append(contract_edges_dif)
        if 'contract_iso' in differentials:
            contract_iso_edges_dif = RestrictedContractEdgesD(
                contract_edges_dif)
            differential_list.append(contract_iso_edges_dif)
            print("Attention: contract_iso operates on nonzero cohomology entries only, so they need to be computed before!")
        super(CHairyGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))

    def print_dim_and_eulerchar(self):
        for h in self.h_range:
            for l in self.l_range:
                ds = [CHairyGraphVS(v, l, h, self.even_edges).get_dimension()
                      for v in self.v_range]
                eul = sum([(1 if j % 2 == 0 else -1) *
                           d for j, d in enumerate(ds)])
                print("Dimensions (h,l) ",
                      h, l, self.sub_type, ":", ds, "Euler", eul)

    def print_cohomology_dim(self):
        for h in self.h_range:
            for l in self.l_range:
                cohomdict = {}
                for v in self.v_range:
                    D1 = ContractEdgesGO.generate_operator(
                        v, l, h, self.even_edges)
                    D2 = ContractEdgesGO.generate_operator(
                        v+1, l, h, self.even_edges)
                    try:
                        d = CHairyGraphVS(
                            v, l, h, self.even_edges).get_dimension()
                        r1 = D1.get_matrix_rank()
                        r2 = D2.get_matrix_rank()
                        cohomdict[v] = d-r1-r2
                    except:
                        pass

                print("Cohomology Dimensions (h,l) ",
                      h, l, self.sub_type, ":", cohomdict)
