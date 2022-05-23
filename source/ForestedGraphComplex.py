"""Graph complexes with a distinguished spanning forest, such as computing the cohomology
of Out(F_n). Hairs are also allowed and distinguishable"""


__all__ = ['graph_type', 'sub_types', 'ForestedGVS', 'ForestedGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
           'ColorEdgesGO', 'ColorEdgesD', 'ForestedGC']

import itertools
from sage.all import *
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import NautyInterface
import Parameters
import SymmetricGraphComplex
import StoreLoad
import Log


logger = Log.logger.getChild('forested_graph_complex')

graph_type = "forested"

sub_types = {True: "even_edges", False: "odd_edges"}


# ------- Graph Vector Space --------
class PreForestedGVS(GraphVectorSpace.GraphVectorSpace):
    """Forested graph vector space.
    Does not implement any symmetry handling, i.e., no graph is zero by symmetry.
    These graph vector spaces are intermediate objects used in the generation algorithm of the
    "actual" forested graph complex ForestedGVS.

    Sub vector space with specified number of vertices, loops, marked edges, hairs
    and and at least trivalent vertices.

    Graphs have a forest of marked edges.
    There may be multiple unmarked edges, but no unmarked tadpoles.
    The marked edges are the ordinary edges between the first n_vertices vertices.
    The unmarked edges are encoded by the following n_unmarked_edges bivalent vertices.
    The following n_hairs vertices correspond to the hairs.
    The edges that form the hairs cannot be marked.

    Attributes:
        - n_vertices (int): Number of vertices.
        - n_loops (int): Number of loops.
        - n_edges (int): Number of edges.

    """

    def __init__(self, n_vertices, n_loops, n_marked_edges, n_hairs):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices, not counting the hairs.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_marked_edges = n_marked_edges
        self.n_hairs = n_hairs
        # n_edges are only the internal edges, not the hairs
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.n_unmarked_edges = self.n_edges - n_marked_edges
        self.sub_type = "pre"
        super(PreForestedGVS, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops \
            and self.n_marked_edges == other.n_marked_edges and self.n_hairs == other.n_hairs

    def __str__(self):
        return ("PreForestedGVS_%s_%s_%s_%s" % self.get_ordered_param_dict().get_value_tuple())

    def __hash__(self):
        return hash(str(self))

    def get_basis_file_path(self):
        s = "gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops),
                                   ('marked_edges', self.n_marked_edges), ('hairs', self.n_hairs)])

    def get_partition(self):
        # All internal vertices are in color 0
        # the unmarked-edge-vertices are in color 1
        # the hair vertices are in colors 2,...,n_hairs+1.
        return [list(range(0, self.n_vertices))] + [list(range(self.n_vertices, self.n_vertices + self.n_unmarked_edges))] \
            + [[j] for j in range(self.n_vertices + self.n_unmarked_edges,
                                  self.n_vertices + self.n_unmarked_edges + self.n_hairs)]

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        # The marked edges can at most form a spanning tree
        return (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs) \
            and self.n_vertices > 0 and self.n_loops >= 0 and self.n_hairs >= 0 \
            and self.n_unmarked_edges >= 0 and self.n_marked_edges >= 0 \
            and self.n_marked_edges <= self.n_vertices - 1
        #          and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2 \

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        # BUT: For algorithmic reasons the graphs with fewer marked edges must be handled first
        if not self.is_valid():
            return 0
        return self.n_vertices ** self.n_hairs \
            * binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices) \
            + self.n_marked_edges

    def get_hairy_graphs(self, nvertices, nloops, nhairs, include_novertgraph=false):
        """ Produces all connected hairy graphs with nhairs hairs, that are the last vertices in the ordering.
        Graphs can have multiple hairs or multiple edges, but not tadpoles.
        Edges are encoded via bivalent vertices of the second color.
        :param include_novertgraph: Whether to include the graph with one edge and no vertices as a two-hair graph
        :type include_novertgraph: bool
        """
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2.
        # Degree 1 vertices are hairs, degree 2 vertices are edges and are removed later.
        nedges = nloops + nvertices - 1   # number of internal edges
        n_vertices_1 = nvertices
        n_vertices_2 = nhairs + nedges    # vertices for edges and hairs
        # number of edges of the bipartite graphs to be generated
        n_edges_bip = nhairs + 2 * nedges
        deg_range_1 = (3, n_edges_bip + 1)
        deg_range_2 = (1, 2)

        # check if valid
        unordered = []
        if (nvertices >= 1 and nloops >= 0 and nhairs >= 0 and n_edges_bip >= n_vertices_2
            and n_edges_bip <= 2*n_vertices_2 and n_edges_bip >= 3 * n_vertices_1
                and n_edges_bip <= n_vertices_1 * n_vertices_2):

            bips = NautyInterface.list_bipartite_graphs3(
                n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip, 2)
            # bips[2].show()
            # print(bips[2].adjacency_matrix())

            # move hair vertices to the end. We are optimistic and assume that the bivalent vertices are always last
            # when canonically ordered
            pp = list(range(nvertices)) + list(range(nvertices+nedges, nvertices +
                                                     nedges+nhairs)) + list(range(nvertices, nvertices+nedges))
            # unordered = (G.relabel(pp, inplace=False) for G in bips)
            for G in bips:
                yield G.relabel(pp, inplace=False)

            # unordered = [self._bip_to_ordinary(
            #     G, nvertices, nedges, nhairs) for G in bipartite_graphs]
        # Produce all permutations of the hairs
        # all_perm = [ range(0,nvertices) + p for p in Permutations(range(nvertices, nvertices+nhairs)) ]
        # return [G.relabel(p, inplace=False) for p in all_perm ]
        if include_novertgraph and nvertices == 0 and nhairs == 2 and nloops == 0:
            # unordered.append(Graph([(0, 1)]))
            yield Graph([(0, 1)])
        # return unordered

    # def get_required_prevs(self):
    #     """Returns the list of PreForestedGVS that need to be built before the present one.
    #     Concretely, these are all the vector spaces with the same parameters but fewer markings. """
    #     for m in range(self.n_marked_edges):
    #         yield PreForestedGVS(self.n_vertices, self.n_loops, m, self.n_hairs)

    def get_generating_graphs(self):
        # Generates all forested graphs.
        # The algorithm is such that if n_marked_edges = 0 one just created all graphs.
        # Otherwise, load the list of graphs with one less marked edges, and mark one.
        if not self.is_valid():
            return

        if self.n_marked_edges == 0:
            yield from self.get_hairy_graphs(self.n_vertices, self.n_loops, self.n_hairs, False)
            return

        VS = PreForestedGVS(self.n_vertices, self.n_loops,
                            self.n_marked_edges-1, self.n_hairs)
        graphs_oneless = VS.get_basis()
        # res = []
        for G in graphs_oneless:
            for i in range(self.n_vertices, self.n_vertices+self.n_unmarked_edges+1):
                nb = G.neighbors(i)
                if len(nb) != 2:
                    raise ValueError(
                        '%s: Vertices of second colour should have 2 neighbours' % str(self))

                if G.has_edge(nb):
                    continue

                GG = copy(G)
                GG.add_edge(nb)
                GG.delete_vertex(i)
                GG.relabel(range(0, GG.order()))

                # check for loop
                if (GG.subgraph(range(self.n_vertices)).is_forest()):
                    yield GG
                    # res.append(GG)

        # return res

    def perm_sign(self, G, p):
        return 1


class ForestedGVS(SymmetricGraphComplex.SymmetricGraphVectorSpace):
    """Forested graph vector space.
    Similar to PreForestedGVS, but implements the symmetry handling.

    Sub vector space with specified number of vertices, loops, marked edges, hairs
    and and at least trivalent vertices.

    Graphs have a forest of marked edges.
    There may be multiple unmarked edges, and unmarked tadpoles.
    The marked edges are the ordinary edges between the first n_vertices vertices.
    The unmarked edges are encoded by the following n_unmarked_edges vertices.
    Bivalent such vertices encode edges, and univalent encode tadpoles.
    The following n_hairs vertices correspond to the hairs.
    The edges that form the hairs cannot be marked.

    We count a tadpole as valence one.
    Hence we miss tadpole vertices of valence three and vertices with two tadpoles.
    This is only problematic if the graph has exactly one vertex, we ignore this case for now.

    Attributes:
        - n_vertices (int): Number of vertices.
        - n_loops (int): Number of loops.
        - n_edges (int): Number of edges.

    """

    def __init__(self, n_vertices, n_loops, n_marked_edges, n_hairs, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: int: Number of vertices, not counting the hairs.
        :type n_vertices: int
        :param n_loops: int: Number of loops.
        :type n_loops: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_marked_edges = n_marked_edges
        self.n_hairs = n_hairs
        # n_edges are only the internal edges, not the hairs
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.n_unmarked_edges = self.n_edges - n_marked_edges
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        self.preVS = PreForestedGVS(
            n_vertices, n_loops, n_marked_edges, n_hairs)

        super(ForestedGVS, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops \
            and self.n_marked_edges == other.n_marked_edges and self.n_hairs == other.n_hairs \
            and self.even_edges == other.even_edges

    def __str__(self):
        return ("ForestedGVS_%s_%s_%s_%s" % self.get_ordered_param_dict().get_value_tuple()) + self.sub_type

    def __hash__(self):
        return hash(str(self))

    def get_basis_file_path(self):
        s = "gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops),
                                   ('marked_edges', self.n_marked_edges), ('hairs', self.n_hairs)])

    def get_partition(self):
        # All internal vertices are in color 0
        # the unmarked-edge-vertices are in color 1
        # the hair vertices are in colors 2,...,n_hairs+1.
        return self.preVS.get_partition()

    def is_valid(self):
        return self.preVS.is_valid()

    def get_work_estimate(self):
        return self.preVS.get_work_estimate()

    def get_generating_graphs(self):
        if not self.is_valid():
            print("Invalid")
            return

        # we assume the basis of the intermediate GVS has alread been constructed
        # We need to add (tp many) tadpoles to graphs and permute hairs

        # res = []
        # we have no tadpoles if edges are even
        maxtp = 0 if self.even_edges else self.n_loops
        for tp in range(0, maxtp+1):
            # newgs = []
            preVS = PreForestedGVS(
                self.n_vertices, self.n_loops-tp, self.n_marked_edges, self.n_hairs+tp)
            preGs = preVS.get_basis()

            if self.n_hairs+tp == 0:
                # newgs = [G for G in preGs]
                yield from preGs
            else:
                # Produce all permutations of the hairs, including those that are encoding tadpoles
                id = Permutation(range(1, tp+1))
                p1s = [p for pp in Permutations(self.n_hairs)
                       for p in id.shifted_shuffle(pp)]
                idv = list(range(0, self.n_vertices+self.n_unmarked_edges-tp))
                all_perm = [idv + [j+self.n_vertices - 1 +
                                   self.n_unmarked_edges-tp for j in p] for p in p1s]

                # all_perm = [list(range(0, self.n_vertices+self.n_unmarked_edges))
                #             + list(p)
                #             for p in itertools.permutations(range(self.n_vertices+self.n_unmarked_edges,
                #                                                   self.n_vertices+self.n_hairs+self.n_unmarked_edges))]

                # newgs = [G.relabel(p, inplace=False)
                #          for G in preGs for p in all_perm]
                for G in preGs:
                    for p in all_perm:
                        yield G.relabel(p, inplace=False)

            # res = res+newgs

        # We count a tadpole as valence one.
        # Hence we miss tadpole vertices of valence three and vertices with two tadpoles.
        # This is only problematic in some edge cases, that we ignore for now.

        # return res

    def get_required_prevs(self):
        """Returns a list of PreForestedGVS that are required to build the basis of this vector space.
        Also includes the PreForestedGVS (...with fewer marked edges) that are needed by those PreForestedGVS
        in turn."""

        if not self.is_valid():
            return

        maxtp = 0 if self.even_edges else self.n_loops
        for tp in range(maxtp+1):
            for m in range(self.n_marked_edges+1):
                yield PreForestedGVS(
                    self.n_vertices, self.n_loops-tp, m, self.n_hairs+tp)

    def label_marked_edges(self, G):
        i = 0
        for (u, v) in G.edges(labels=False):
            if (u < self.n_vertices and v < self.n_vertices):
                G.set_edge_label(u, v, i)
                i += 1

    def perm_sign(self, G, p):
        if self.even_edges:
            # The total sign is:
            # (a:induced sign on vertices)
            # * (b:induced sign of the unmarked-edge-permutation)
            # * (c:induced sign unmarked-edge orientations)
            # * (d:induced sign marked-edge orientations)
            # note that we have no tadpoles in case of even edges

            # This is a*b
            sign = Shared.Perm(
                [j for j in p if (j < self.n_vertices+self.n_unmarked_edges)]).signature()
            # Now  * d
            for (u, v) in G.edges(labels=False):
                if (u < self.n_vertices and v < self.n_vertices):
                    # We assume the edge is always directed from the larger to smaller index
                    if (u < v and p[u] > p[v]) or (u > v and p[u] < p[v]):
                        sign *= -1
            # Finally * c
            for i in range(self.n_vertices, self.n_vertices+self.n_unmarked_edges):
                nb = G.neighbors(i)
                if len(nb) != 2:
                    # This is a graph with a tadpole, and hence zero... we return zero for now, although this is not ideal.
                    return 0
                    # raise ValueError(
                    #    '%s: Vertices of second colour should have 2 neighbours' % str(self))
                u = nb[0]
                v = nb[1]
                if (u < v and p[u] > p[v]) or (u > v and p[u] < p[v]):
                    sign *= -1
            return sign
        else:
            # The sign is (induced sign of the marked-edge-permutation)
            # We assume the edges are always lexicographically ordered
            # For the computation we use that G.edges() returns the edges in lex ordering
            # We first label the edges on a copy of G lexicographically
            G1 = copy(G)
            self.label_marked_edges(G1)
            # print("edges before ", [j for (u, v, j) in G1.edges()])
            # We permute the graph, and read of the new labels, ...
            # but only those between internal vertices of the first color, i.e., the first n_vertices ones
            G1.relabel(p, inplace=True)
            # print("edges after ", [(u, v, j) for (u, v, j) in G1.edges()])

            return Shared.Perm([j for (u, v, j) in G1.edges() if (u < self.n_vertices and v < self.n_vertices)]).signature()

    def get_n(self):
        return self.n_hairs

    def vertex_permutation_from_permutation(self, p):
        return list(range(0, self.n_vertices)) + [j+self.n_vertices-1 for j in p]

    def get_isotypical_projector(self, rep_index):
        return SymmProjector(self, rep_index)


class ForestedGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of forested graph vector spaces with specified edge parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, m_range, h_range, even_edges):
        """Initialize the sum vector space.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param m_range: Range for the number of marked edges.
        :type m_range: range
        :param h_range: Range for the number of hairs.
        :type h_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.v_range = v_range
        self.l_range = l_range
        self.m_range = m_range
        self.h_range = h_range

        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        vs_list = [ForestedGVS(v, l, m, h, self.even_edges) for (
            v, l, m, h) in itertools.product(self.v_range, self.l_range, self.m_range, self.h_range)]

        super(ForestedGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('marked_edges', self.m_range), ('hairs', self.h_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s_%s" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


class PreForestedGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of forested graph vector spaces with specified edge parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, m_range, h_range):
        """Initialize the sum vector space.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param m_range: Range for the number of marked edges.
        :type m_range: range
        :param h_range: Range for the number of hairs.
        :type h_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.v_range = v_range
        self.l_range = l_range
        self.m_range = m_range
        self.h_range = h_range

        self.sub_type = "pre"

        vs_list = [PreForestedGVS(v, l, m, h) for (
            v, l, m, h) in itertools.product(self.v_range, self.l_range, self.m_range, self.h_range)]

        super(PreForestedGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s pre graphs' % (graph_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('marked_edges', self.m_range), ('hairs', self.h_range)])

    def get_info_plot_path(self):
        s = "info_pre_vector_space_%s_%s" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


class PreForestedGraphSumVS2(GraphVectorSpace.SumVectorSpace):
    """This is for holding an arbitrary list of PreForestedGVS that are not necessarily in a consecutive range
     in the parameter table."""

    def __init__(self, vs_list):
        self.sub_type = "pre"
        super(PreForestedGraphSumVS2, self).__init__(vs_list)

    def get_type(self):
        return f'{graph_type} pre graphs'

    def get_ordered_param_range_dict(self):
        return {"parameters": "various"}

    def get_info_plot_path(self):
        s = "info_pre_vector_space_%s_%s" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


# ------- Operators --------
class ContractEdgesGO(SymmetricGraphComplex.SymmetricGraphOperator):
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
        super(ContractEdgesGO, self).__init__(domain, target)

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
            and domain.n_marked_edges - 1 == target.n_marked_edges \
            and domain.n_hairs == target.n_hairs \
            and domain.even_edges == target.even_edges

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_marked_edges, n_hairs, even_edges):
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
        domain = ForestedGVS(n_vertices, n_loops,
                             n_marked_edges, n_hairs, even_edges)
        target = ForestedGVS(n_vertices - 1, n_loops,
                             n_marked_edges-1, n_hairs, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        # Operates on the graph G by contracting a marked edge and unifying the adjacent vertices.
        image = []
        for (u, v) in G.edges(labels=False):
            # only contract marked edges
            if (u < self.domain.n_vertices and v < self.domain.n_vertices):
                # move the two vertices to be merged to the first two positions
                pp = Shared.permute_to_left(
                    (u, v), range(0, G.order()))
                sgn = self.domain.perm_sign(G, pp)
                G1 = copy(G)
                G1.relabel(pp, inplace=True)
                self.domain.label_marked_edges(G1)
                # previous_size = G1.size()
                G1.merge_vertices([0, 1])
                # if (previous_size - G1.size()) != 1:
                #     continue
                G1.relabel(list(range(0, G1.order())), inplace=True)
                if not self.domain.even_edges:
                    # for odd edges compute the sign of the permutation of internal edges
                    p = [j for (a, b, j) in G1.edges() if (
                        a < self.target.n_vertices and b < self.target.n_vertices)]
                    # If we removed more then one marked edge stop
                    if len(p) < self.domain.n_marked_edges-1:
                        print("This should not happen....")
                        continue
                    sgn *= Permutation(p).signature()
                else:
                    # There is no further sign for even edges
                    sgn *= 1  # TODO overall sign for even edges
                image.append((G1, sgn))
        return image

    def restrict_to_isotypical_component(self, rep_index):
        return RestrictedContractEdgesGO(self, rep_index)


class RestrictedContractEdgesGO(SymmetricGraphComplex.SymmetricRestrictedOperatorMatrix):

    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d_r%d.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_%d_r%d_rank.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def is_match(self, domain, target):
        return ContractEdgesGO.is_match(domain.vs, target.vs) and domain.rep_index == target.rep_index


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: OrdinaryGraphSumVS
        """
        super(ContractEdgesD, self).__init__(sum_vector_space,
                                             ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_D_%s_%s" % (graph_type, sub_type)
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

    def get_ordered_param_dict2(self):
        do = self.domain
        return Shared.OrderedDict([('vertices', do.n_vertices), ('loops', do.n_loops), ('marked edges', do.n_marked_edges), ('hairs', do.n_hairs), ('rep_index', self.rep_index)])

    def get_matrix_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d_rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "projectionO%d_%d_%d_%d_%d.txt.rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)


# class UnmarkEdgesD(GraphOperator.Differential):
#     """Unmark edges differential. Sums over marked edges and unmarks one."""

#     def __init__(self, sum_vector_space):
#         """Initialize the contract edges differential with the underlying sum vector space.

#         :param sum_vector_space: Underlying vector space.
#         :type sum_vector_space: OrdinaryGraphSumVS
#         """
#         super(UnmarkEdgesD, self).__init__(sum_vector_space,
#                                            ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

#     def get_type(self):
#         return 'contract edges'

#     def get_cohomology_plot_path(self):
#         sub_type = self.sum_vector_space.sub_type
#         s = "cohomology_dim_contrct_edges_D_%s_%s" % (graph_type, sub_type)
#         return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

#     def get_info_plot_path(self):
#         sub_type = self.sum_vector_space.sub_type
#         s = "info_contract_edges_D_%s_%s" % (graph_type, sub_type)
#         return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class UnmarkEdgesGO(SymmetricGraphComplex.SymmetricGraphOperator):
    """Unmark edges graph operator.
    Makes one marked edge into an unmarked edge.

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
        if not UnmarkEdgesGO.is_match(domain, target):
            raise ValueError(
                "Domain and target not consistent for delete edges operator")
        self.sub_type = domain.sub_type
        super(UnmarkEdgesGO, self).__init__(domain, target)

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
        return domain.n_vertices == target.n_vertices and domain.n_loops == target.n_loops \
            and domain.n_hairs == target.n_hairs and domain.n_marked_edges-1 == target.n_marked_edges \
            and domain.even_edges == target.even_edges

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_marked_edges, n_hairs, even_edges):
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
        domain = ForestedGVS(n_vertices, n_loops,
                             n_marked_edges, n_hairs, even_edges)
        target = ForestedGVS(n_vertices, n_loops,
                             n_marked_edges-1, n_hairs, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "unmark_edgesD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "unmark_edgesD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "unmark_edgesD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "unmark_edgesD%d_%d_%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        return 'unmark edges'

    def operate_on(self, G):
        # Operates on the graph G by unmarking a marked edge.
        image = []
        # create a new graph with one more vertex of second color that encodes the new unmarked edge
        # the new vertex will be the first among the second color vertices (index n_vertices)
        GG = copy(G)
        GG.add_vertex()
        GG.relabel(list(range(0, self.domain.n_vertices)) + list(range(self.domain.n_vertices +
                   1, GG.order())) + [self.domain.n_vertices], inplace=True)
        i = 0  # counts the index of the edge to be unmarked for sign purposes
        for (u, v) in G.edges(labels=False):
            # only unmark marked edges...
            if (u < self.domain.n_vertices and v < self.domain.n_vertices):
                G1 = copy(GG)
                G1.delete_edge((u, v))
                G1.add_edge((u, self.domain.n_vertices))
                G1.add_edge((v, self.domain.n_vertices))
                i += 1
                # For even edges the sign is given by nr 0f vertices
                sgn = 1
                if self.domain.even_edges:
                    sgn = 1 if (self.domain.n_vertices % 2 == 0) else -1
                if not self.domain.even_edges:
                    sgn = 1 if (i % 2 == 0) else -1
                image.append((G1, sgn))
        return image

    def restrict_to_isotypical_component(self, rep_index):
        return RestrictedUnmarkEdgesGO(self, rep_index)


class RestrictedUnmarkEdgesGO(SymmetricGraphComplex.SymmetricRestrictedOperatorMatrix):

    def get_matrix_file_path(self):
        s = "unmark_edgesD%d_%d_%d_%d_r%d.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_rank_file_path(self):
        s = "unmark_edgesD%d_%d_%d_%d_r%d_rank.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def is_match(self, domain, target):
        return UnmarkEdgesGO.is_match(domain.vs, target.vs) and domain.rep_index == target.rep_index


class UnmarkEdgesD(GraphOperator.Differential):
    """Unmark edges differential.
    """

    def __init__(self, sum_vector_space):
        """Initialize the delete edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: OrdinaryGraphSumVS
        """
        super(UnmarkEdgesD, self).__init__(sum_vector_space,
                                           UnmarkEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'unmark edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_unmark_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_unmark_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class RestrictedUnmarkEdgesD(SymmetricGraphComplex.SymmetricDifferential):
    def get_type(self):
        return 'isotypical unmark edges'

    def get_cohomology_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "cohomology_dim_unmark_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "info_unmark_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


# ------- Graph Complexes --------
class ForestedGC(GraphComplex.GraphComplex):
    """Graph complex for forested graphs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, m_range, h_range, even_edges, differentials):
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
        self.m_range = m_range
        self.h_range = h_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        sum_vector_space = ForestedGraphSumVS(
            v_range, l_range, m_range, h_range, even_edges)
        differential_list = []
        if not set(differentials) <= {'contract', 'unmark'}:
            raise ValueError(
                "Differentials for forested graph complex: 'contract', 'unmark'")
        if 'contract' in differentials or 'contractiso' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            if 'contract' in differentials:
                differential_list.append(contract_edges_dif)
            if 'contract_iso' in differentials:
                differential_list.append(
                    RestrictedContractEdgesD(contract_edges_dif))
        if 'unmark' in differentials or 'unmark_iso' in differentials:
            delete_edges_dif = UnmarkEdgesD(sum_vector_space)
            if 'unmark' in differentials:
                differential_list.append(delete_edges_dif)
            if 'unmark_iso' in differentials:
                differential_list.append(
                    RestrictedUnmarkEdgesD(delete_edges_dif))
        super(ForestedGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))

    def compute_all_pregraphs(self, **kwargs):
        print("Determining and building required pre-vs:")
        vsset = {
            prevs for vs in self.sum_vector_space.vs_list for prevs in vs.get_required_prevs()}
        for vs in vsset:
            print(vs)

        sumvs = PreForestedGraphSumVS2(list(vsset))
        sumvs.build_basis(**kwargs)

    def build_basis(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        print("Building auxiliary pregraphs...")
        self.compute_all_pregraphs(ignore_existing_files=ignore_existing_files,
                                   n_jobs=n_jobs, progress_bar=progress_bar, info_tracker=info_tracker)
        # PreForestedGraphSumVS.compute_all_pregraphs(max(self.v_range),
        # max(self.l_range), max(self.m_range), max(self.h_range), self.even_edges, ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar, info_tracker=info_tracker)
        print("Done.")
        return super().build_basis(ignore_existing_files, n_jobs, progress_bar, info_tracker)


# ------------- Bicomplex ------------------------

class ContractUnmarkBiOM(SymmetricGraphComplex.SymmetricBiOperatorMatrix):
    """Bi operator matrix based on the differentials contract edges and unmark edges.

    Attributes:
            - sub_type (str): Sub type of graphs.
    """

    def __init__(self, domain, target):
        self.sub_type = domain.sub_type
        super(ContractUnmarkBiOM, self).__init__(domain, target, ContractEdgesGO,
                                                 UnmarkEdgesGO)

    @classmethod
    def generate_operator(cls, n_loops, n_marked_edges, n_hairs, even_edges):
        domain = ForestedDegSlice(n_loops,
                                  n_marked_edges, n_hairs, even_edges)
        target = ForestedDegSlice(n_loops,
                                  n_marked_edges-1, n_hairs, even_edges)
        return cls(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target degree slices match to generate a corresponding bi operator matrix.

        The bi operator reduces the degree by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: ForestedDegSlice
        :param target: Potential target vector space of the operator.
        :type target: ForestedDegSlice
        :return: bool: True if domain and target match to generate a corresponding bi operator matrix.
        :rtype: bool
        """
        return domain.n_marked_edges - 1 == target.n_marked_edges \
            and domain.n_loops == target.n_loops \
            and domain.n_hairs == target.n_hairs \
            and domain.even_edges == target.even_edges

    def get_matrix_file_path(self):
        s = "bi_D_contract_unmark_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_contract_unmark_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def restrict_to_isotypical_component(self, rep_index):
        return RestrictedContractUnmarkGO(self, rep_index)


class RestrictedContractUnmarkGO(SymmetricGraphComplex.SymmetricRestrictedOperatorMatrix):

    def get_matrix_file_path(self):
        s = "bi_D_contract_unmark_%d_%d_%d_r%d.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_contract_unmark_%d_%d_%d_r%d_rank.txt" % (
            self.domain.vs.get_ordered_param_dict().get_value_tuple() + (self.rep_index,))
        return os.path.join(Parameters.data_dir, graph_type, self.opD.sub_type, s)

    def is_match(self, domain, target):
        return ContractUnmarkBiOM.is_match(domain.vs, target.vs) and domain.rep_index == target.rep_index


class ForestedDegSlice(SymmetricGraphComplex.SymmetricDegSlice):
    """Degree slice of forested graphs

    Total degree = marked edges

    Attributes:
        - even_edges (bool): True for even edges, False for odd edges.

    """

    def is_complete(self):
        for vs in self.vs_list:
            if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
                return False
        return True

    def __init__(self, n_loops, n_marked_edges, n_hairs, even_edges):
        """Initialize the degree slice.

        :param deg: Total degree of the degree slice.
        :type deg: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_loops = n_loops
        self.n_marked_edges = n_marked_edges
        self.n_hairs = n_hairs
        self.even_edges = even_edges
        self.sub_type = sub_types.get(even_edges)
        max_vertices = 2*n_loops-2 + n_hairs
        min_vertices = n_marked_edges+1
        super(ForestedDegSlice, self).__init__(
            [ForestedGVS(v, n_loops, n_marked_edges, n_hairs, even_edges)
             for v in range(min_vertices, max_vertices + 1)],
            n_marked_edges)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('loops', self.n_loops), ('marked_edges', self.n_marked_edges), ('hairs', self.n_hairs)])

    def __eq__(self, other):
        return self.n_loops == other.n_loops \
            and self.n_marked_edges == other.n_marked_edges and self.n_hairs == other.n_hairs

    def __str__(self):
        return ("ForestedDegSlice_%s_%s_%s" % self.get_ordered_param_dict().get_value_tuple()) + self.sub_type

    def __hash__(self):
        return hash(str(self))

    def get_info_plot_path(self):
        s = "info_vertex_loop_degree_slice_deg_%d_%d_%d_%s_%s" % (
            self.n_loops, self.n_marked_edges, self.n_hairs, graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

    def get_n(self):
        return self.n_hairs

    def get_isotypical_projector(self, rep_index):
        """Returns the SymmetricProjectionOperator corresponding to the isotypical component corresponding to
        the rep_index-th irrep (as in Partitions(n))
        """
        return SymmProjectorDegSlice(self, rep_index)


class SymmProjectorDegSlice(SymmetricGraphComplex.SymmetricProjectionOperatorDegSlice):
    def __init__(self, domain, rep_index):
        """
        : param domain: Domain vector space of the operator.
        : type domain: HairyGraphVS
        : param rep_index: The index of the representation in the list produced by Partitions(h).
        : type rep_index: int
        """
        self.sub_type = domain.sub_type

        super(SymmProjectorDegSlice, self).__init__(domain, rep_index)

    def get_ordered_param_dict2(self):
        do = self.domain
        return Shared.OrderedDict([('loops', do.n_loops), ('marked edges', do.n_marked_edges), ('hairs', do.n_hairs), ('rep_index', self.rep_index)])

    def get_matrix_file_path(self):
        s = "projectionODegSlice%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "projectionODegSlice%d_%d_%d_%d_rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "projectionODegSlice%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "projectionODegSlice%d_%d_%d_%d.txt.rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)


class ForestedBigradedSumVS(GraphVectorSpace.SumVectorSpace):
    """Bi graded vector space based on ordinary simple graphs.

    Bi grading according to the number of vertices and loops.
    Direct sum of degree slices.

    Attributes:
        - deg_range (range): Range for the total degree.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, l_range, m_range, h_range, even_edges):
        """ Initialize the bi graded vector space.

        :param deg_range: Range for the degree.
        :type deg_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.l_range = l_range
        self.m_range = m_range
        self.h_range = h_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(even_edges)
        super(ForestedBigradedSumVS, self).__init__([ForestedDegSlice(l, m, h, self.even_edges)
                                                     for l in l_range for m in m_range for h in h_range])

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('loops', self.l_range), ('marked_edges', self.m_range), ('hairs', self.h_range)])

    def get_info_plot_path(self):
        s = "info_bigraded_vector_space_%s_%s" % (
            graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


class ContractUnmarkD(GraphOperator.Differential):
    """Differential on the bi graded vector space based on the operators contract edges and delete edges.

    Only for graphs with odd edges.
    """

    def __init__(self, graded_sum_vs):
        """Initialize the contract and delete edges differential with the underlying bi graded vector space.

        :param graded_sum_vs: Underlying bi graded vector space.
        :type graded_sum_vs: VertexLoopBigradedSumVS
        """
        super(ContractUnmarkD, self).__init__(graded_sum_vs,
                                              ContractUnmarkBiOM.generate_op_matrix_list(graded_sum_vs))

    def get_type(self):
        return 'contract edges and unmark edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_unmark_edges_D_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_unmark_edges_D_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        s = self.sum_vector_space
        return Shared.OrderedDict([('loops', s.l_range), ('marked_edges', s.m_range), ('hairs', s.h_range)])


class RestrictedContractUnmarkD(SymmetricGraphComplex.SymmetricDifferential):
    def get_type(self):
        return 'isotypical contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_unmark_edges_D_iso_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "info_contract_edges_unmark_edges_D_iso_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class ForestedContractUnmarkBiGC(GraphComplex.GraphComplex):
    """Bi complex based on ordinary simple graphs and the differentials contract edges and unmark edges.

    Attributes:
        - deg_range (range): Range for the total degree.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, l_range, m_range, h_range, even_edges, isotypical=False):
        """Initialize the bi complex.

        :param deg_range: Range for the degree.
        :type deg_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.l_range = l_range
        self.m_range = m_range
        self.h_range = h_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)
        graded_sum_vs = ForestedBigradedSumVS(
            l_range, m_range, h_range, self.even_edges)
        self.contract_unmarkD = ContractUnmarkD(graded_sum_vs)
        if isotypical:
            super(ForestedContractUnmarkBiGC, self).__init__(
                graded_sum_vs, [RestrictedContractUnmarkD(self.contract_unmarkD)])
        else:
            super(ForestedContractUnmarkBiGC, self).__init__(
                graded_sum_vs, [self.contract_unmarkD])

    def __str__(self):
        return '<%s graphs bi-complex with %s>' % (graph_type, str(self.sub_type))

    def print_dim_and_eulerchar(self):
        for h in self.h_range:
            for l in self.l_range:
                ds = [ForestedDegSlice(l, m, h, self.even_edges).get_dimension()
                      for m in self.m_range]
                eul = sum([(1 if j % 2 == 0 else -1) *
                           d for j, d in enumerate(ds)])
                print("Dimensions (h,l) ",
                      h, l, self.sub_type, ":", ds, "Euler", eul)

    def print_cohomology_dim(self):
        for h in self.h_range:
            for l in self.l_range:
                cohomdict = {}
                for m in self.m_range:
                    D1 = ContractUnmarkBiOM.generate_operator(
                        l, m, h, self.even_edges)
                    D2 = ContractUnmarkBiOM.generate_operator(
                        l, m+1, h, self.even_edges)
                    try:
                        d = ForestedDegSlice(
                            l, m, h, self.even_edges).get_dimension()
                        r1 = D1.get_matrix_rank()
                        r2 = D2.get_matrix_rank()
                        cohomdict[m] = d-r1-r2
                    except:
                        pass

                print("Cohomology Dimensions (h,l) ",
                      h, l, self.sub_type, ":", cohomdict)

    def build_basis(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        # print("Building auxiliary pregraphs...")
        # PreForestedGraphSumVS.compute_all_pregraphs(-1,
        #                                             max(self.l_range), max(self.m_range), max(self.h_range), self.even_edges, ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar, info_tracker=info_tracker)
        # print("Done.")
        return super().build_basis(ignore_existing_files, n_jobs, progress_bar, info_tracker)


###### Top vertex number parts ###

class ForestedTopDegSlice(SymmetricGraphComplex.SymmetricDegSlice):
    """Represents the top and top-1-piece of a forested deg slice, i.e.,
    just the graphs with all trivalent, or with one 4-valent vertex.
    """

    def is_complete(self):
        for vs in self.vs_list:
            if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
                return False
        return True

    def __init__(self, n_loops, n_marked_edges, n_hairs, even_edges, top_n):
        """Initialize the degree slice.
        There will be top_n gvs in the deg slice. (only 1 and 2 are used)

        :param deg: Total degree of the degree slice.
        :type deg: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.n_loops = n_loops
        self.n_marked_edges = n_marked_edges
        self.n_hairs = n_hairs
        self.even_edges = even_edges
        self.sub_type = sub_types.get(even_edges)
        self.top_n = top_n
        max_vertices = 2*n_loops-2 + n_hairs
        self.n_vertices = max_vertices
        min_vertices = max_vertices-top_n+1
        super(ForestedTopDegSlice, self).__init__(
            [ForestedGVS(v, n_loops, n_marked_edges, n_hairs, even_edges)
             for v in range(min_vertices, max_vertices + 1)],
            n_marked_edges)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('loops', self.n_loops), ('marked_edges', self.n_marked_edges), ('hairs', self.n_hairs), ('topn', self.top_n)])

    def __eq__(self, other):
        return self.n_loops == other.n_loops \
            and self.n_marked_edges == other.n_marked_edges and self.n_hairs == other.n_hairs \
            and self.even_edges == other.even_edges

    def __str__(self):
        return ("ForestedTopDegSlice_%s_%s_%s_%s" % self.get_ordered_param_dict().get_value_tuple()) + self.sub_type

    def __hash__(self):
        return hash(str(self))

    def get_info_plot_path(self):
        s = "info_vertex_loop_top_degree_slice_deg_%d_%d_%d_%d_%s_%s" % (
            self.n_loops, self.n_marked_edges, self.n_hairs, self.top_n, graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

    def get_n(self):
        return self.n_hairs

    def get_isotypical_projector(self, rep_index):
        """Returns the SymmetricProjectionOperator corresponding to the isotypical component corresponding to
        the rep_index-th irrep (as in Partitions(n))
        """
        return SymmProjectorDegSlice(self, rep_index)


class ContractUnmarkTopBiOM(SymmetricGraphComplex.SymmetricBiOperatorMatrix):
    """
    Represents the stacked matrix
    M=[ contract; unmark ]
    that can be used to compute the rank of unmark U, restricted to the kernel of contract C, as
    rank(U|_{Ker C} ) = rank(M) - rank(C)

    This is used only for the top number of vertices components (trivalent graphs) in each loop order.
    """

    def __init__(self, domain, target):
        self.sub_type = domain.sub_type
        super(ContractUnmarkTopBiOM, self).__init__(domain, target, ContractEdgesGO,
                                                    UnmarkEdgesGO)

    @classmethod
    def generate_operator(cls, n_loops, n_marked_edges, n_hairs, even_edges):
        domain = ForestedTopDegSlice(n_loops,
                                     n_marked_edges, n_hairs, even_edges, 1)
        target = ForestedTopDegSlice(n_loops,
                                     n_marked_edges-1, n_hairs, even_edges, 2)
        return cls(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target degree slices match to generate a corresponding bi operator matrix.

        The bi operator reduces the degree by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: ForestedDegSlice
        :param target: Potential target vector space of the operator.
        :type target: ForestedDegSlice
        :return: bool: True if domain and target match to generate a corresponding bi operator matrix.
        :rtype: bool
        """
        return domain.n_marked_edges - 1 == target.n_marked_edges \
            and domain.n_loops == target.n_loops \
            and domain.n_hairs == target.n_hairs \
            and domain.even_edges == target.even_edges

    def get_matrix_file_path(self):
        s = "bi_D_contract_unmark_top_%d_%d_%d.txt" % self.domain.get_ordered_param_dict(
        ).get_value_tuple()[:-1]
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_contract_unmark_top_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict(
        ).get_value_tuple()[:-1]
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_true_rank(self):
        """ Returns the rank of the restriction of the unmark operator to the kernel of the
        contraction operator.
        """
        Dc = ContractEdgesGO.generate_operator(
            self.domain.n_vertices, self.domain.n_loops, self.domain.n_marked_edges, self.domain.n_hairs, self.domain.even_edges)

        return self.get_matrix_rank() - Dc.get_matrix_rank()

    def restrict_to_isotypical_component(self, rep_index):
        pass
        # return RestrictedContractUnmarkGO(self, rep_index)


class ForestedGraphTopSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of forested graph vector spaces with specified edge parity.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - even_edges (bool): True for even edges, False for odd edges.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, l_range, m_range, h_range, even_edges):
        """Initialize the sum vector space.

        :param v_range: Range for the number of vertices.
        :type v_range: range
        :param l_range: Range for the number of loops.
        :type l_range: range
        :param m_range: Range for the number of marked edges.
        :type m_range: range
        :param h_range: Range for the number of hairs.
        :type h_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        """
        self.l_range = l_range
        self.m_range = m_range
        self.h_range = h_range

        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        vs_list = [ForestedGVS(2*l-2+h, l, m, h, self.even_edges) for (
            l, m, h) in itertools.product(self.l_range, self.m_range, self.h_range)] \
            + \
            [ForestedGVS(2*l-2+h-1, l, m, h, self.even_edges) for (
                l, m, h) in itertools.product(self.l_range, self.m_range, self.h_range)]

        super(ForestedGraphTopSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('loops', self.l_range), ('marked_edges', self.m_range), ('hairs', self.h_range)])

    def compute_all_pregraphs(self, **kwargs):
        print("Determining and building required pre-vs:")
        vsset = {
            prevs for vs in self.vs_list for prevs in vs.get_required_prevs()}
        for vs in vsset:
            print(vs)

        sumvs = PreForestedGraphSumVS2(list(vsset))
        sumvs.build_basis(**kwargs)

    def get_info_plot_path(self):
        s = "info_vector_space_top_%s_%s" % (graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


class ContractUnmarkTopD(GraphOperator.Differential):
    """ Represents a collection of ContractUnmarkTopBiOM.
    This class is also used to compute cohomology.
    """

    def __init__(self, l_range, m_range, h_range, even_edges):
        self.l_range = l_range
        self.m_range = m_range
        self.h_range = h_range
        self.even_edges = even_edges
        op_list = [ContractUnmarkTopBiOM.generate_operator(l, m, h, even_edges)
                   for l in l_range
                   for m in m_range
                   for h in h_range]
        sum_vs = ForestedGraphTopSumVS(l_range, m_range, h_range, even_edges)
        super(ContractUnmarkTopD, self).__init__(sum_vs, op_list)

    def get_type(self):
        return 'contract edges and unmark edges top'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_unmark_edges_top_D_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_unmark_edges_top_D_%s_%s" % (
            graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        s = self.sum_vector_space
        return Shared.OrderedDict([('loops', s.l_range), ('marked_edges', s.m_range), ('hairs', s.h_range)])

    def _get_single_cohomology(self, n_loops, n_marked, n_hairs):
        """
        Computes a single cohomology dimension.
        """
        opD = ContractUnmarkTopBiOM.generate_operator(
            n_loops, n_marked, n_hairs, self.even_edges)
        opDD = ContractUnmarkTopBiOM.generate_operator(
            n_loops, n_marked+1, n_hairs, self.even_edges)
        n_vertices = opD.domain.n_vertices

        opC = ContractEdgesGO.generate_operator(
            n_vertices, n_loops, n_marked, n_hairs, self.even_edges)

        try:
            # dimension of the kernel of the contraction
            dimV = opD.get_domain().get_dimension() - opC.get_matrix_rank()
        except StoreLoad.FileNotFoundError:
            logger.info(
                "Cannot compute cohomology: First build basis for %s " % str(opD.get_domain()))
            return None
        if dimV == 0:
            return '*'
        if opD.is_valid():
            try:
                rankD = opD.get_true_rank()
            except StoreLoad.FileNotFoundError:
                logger.info(
                    "Cannot compute cohomology: Matrix rank not calculated for %s " % str(opD))
                return None
        else:
            rankD = 0
        if opDD.is_valid():
            try:
                rankDD = opDD.get_true_rank()
            except StoreLoad.FileNotFoundError:
                logger.info(
                    "Cannot compute cohomology: Matrix rank not calculated for %s " % str(opDD))
                return None
        else:
            rankDD = 0
        cohomology_dim = dimV - rankD - rankDD
        if cohomology_dim < 0:
            raise ValueError("Negative cohomology dimension for %s (%d - %d - %d)" %
                             (str(opD.domain), dimV, rankD, rankDD))
            # logger.error("Negative cohomology dimension for %s" % str(opD.domain))
        return cohomology_dim

    def get_cohomology_dim_dict(self):
        ms = list(self.m_range)[:-1]
        return {(l, m, h): self._get_single_cohomology(l, m, h)
                for l in self.l_range
                for h in self.h_range
                for m in ms}

        # return super().get_cohomology_dim_dict()

    def build_basis(self, **kwargs):
        self.sum_vector_space.compute_all_pregraphs(**kwargs)
        self.sum_vector_space.build_basis(**kwargs)

    def compute_rank(self, sage=None, linbox=None, rheinfall=None, sort_key='size', ignore_existing_files=False, n_jobs=1, info_tracker=False):
        # compute ranks of contractto operators
        super().compute_rank(sage, linbox, rheinfall, sort_key,
                             ignore_existing_files, n_jobs, info_tracker)
        # compute ranks of contract operators that are also necessary to have
        print("Computing contract operator ranks...")
        coplist = [ContractEdgesGO.generate_operator(2*l-2+h, l, m, h, self.even_edges)
                   for l in self.l_range
                   for m in self.m_range
                   for h in self.h_range]
        oc = GraphOperator.OperatorMatrixCollection(
            self.sum_vector_space, coplist)
        oc.compute_rank(sage, linbox, rheinfall, sort_key,
                        ignore_existing_files, n_jobs, info_tracker)
