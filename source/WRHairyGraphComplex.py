"""Graph complexes based on simple graphs with numbered hairs and hairs of two decorations.
Implemented Differentials: Contract edges.

This is a reduced version of the WHairy graph complex, defined by setting graphs to zero that have connected components 
with at least one vertex and no omega-hair.

Graphs are realized as simple graphs with 2 extra vertices for epsilon and omega, (index self.n_vertices and self.n_vertices+1).
WARNING: If there is a tadpole the corresponding loop is not part of the graph--one can determine the presence of the tadpole 
from the overall one too small loop number.
TODO: Take care that this does not produce problems
"""


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

graph_type = "wrhairy"

zero_hairs = False      # Option to include zero hairs in the hairy graph complexes.

def dict_to_list(d, n):
    return [(d[j] if j in d else j)   for j in range(n)]
    

# ------- Graph Vector Space --------
class WRHairyGraphVS(GraphVectorSpace.GraphVectorSpace):
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

    def __init__(self, n_vertices, n_loops, n_hairs, n_ws):
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
        self.n_ws = n_ws
        # we count only the internal edges and omega and eps edges, but not the hair edges
        self.n_edges = self.n_loops + self.n_vertices
        self.sub_type = "wr"
        super(WRHairyGraphVS, self).__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(self.n_vertices + self.n_hairs+2, self.n_loops, False)

    def get_type(self):
        return 'wrgraphs'

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs and self.n_ws == other.n_ws

    def __hash__(self):
        return hash("wrgra%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wrgra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "wrgra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('hairs', self.n_hairs), ('ws', self.n_ws)])

    def get_partition(self):
        # All internal vertices are in color 1, the single eps vertex in color 2, the w vertex in color 3
        # and the hair vertices are in colors 4,...,n+3.
        return [list(range(0, self.n_vertices)) ] + [ [j] for j in range(self.n_vertices, self.n_vertices + self.n_hairs+2)]

    def plot_graph(self, G):
        GG = Graph(G, loops=True)
        # add proper tadpole if needed
        if GG.size() < self.n_edges + self.n_hairs:
            GG.add_edge(self.n_vertices, self.n_vertices)
        
        return GG.plot(partition=self.get_partition(), vertex_labels=True)

    def is_valid(self):
        # At least trivalent internal vertices.
        l = (3 * self.n_vertices + self.n_ws <= 2 * self.n_edges + self.n_hairs )
        # Nonnegative number of vertices, non negative number of loops, non-negative or positive number of hairs, 
        # and at least one omega hair.
        l = l and self.n_vertices >= 0 and self.n_loops >= 0 and self.n_hairs >= 0 and self.n_ws >= 1
        # At most a full graph.
        l = l and self.n_edges <= (self.n_vertices+2) * (self.n_vertices + 1) / 2
        return l

    def get_work_estimate(self):
        # TODO
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)



    def get_hairy_graphs(self, nvertices, nloops, nhairs, include_novertgraph=false):
        """ Produces all connected hairy graphs with nhairs hairs, that are the last vertices in the ordering.
        Graphs can have multiple hairs, but not tadpoles or multiple edges.
        :param include_novertgraph: Whether to include the graph with one edge and no vertices as a two-hair graph
        :type include_novertgraph: bool
        """
        # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2.
        # Degree 1 vertices are hairs, degree 2 vertices are edges and are removed later.
        nedges = nloops + nvertices - 1 # number of internal edges
        n_vertices_1 = nvertices
        n_vertices_2 = nhairs + nedges
        n_edges_bip = nhairs + 2 * nedges
        deg_range_1 = (3, n_edges_bip + 1)
        deg_range_2 = (1, 2)

        # check if valid
        unordered=[]
        if (nvertices >= 1 and nloops >= 0 and nhairs >= 0 and n_edges_bip >= n_vertices_2 
                           and n_edges_bip <= 2*n_vertices_2 and n_edges_bip >= 3* n_vertices_1
                           and n_edges_bip <= n_vertices_1 * n_vertices_2):
            bipartite_graphs = NautyInterface.list_bipartite_graphs2(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges_bip)
            unordered = [self._bip_to_ordinary(G,nvertices,nedges,nhairs) for G in bipartite_graphs]
        # Produce all permutations of the hairs
        #all_perm = [ range(0,nvertices) + p for p in Permutations(range(nvertices, nvertices+nhairs)) ]
        #return [G.relabel(p, inplace=False) for p in all_perm ]
        if include_novertgraph and nvertices == 0 and nhairs==2 and nloops == 0:
            unordered.append(Graph([(0,1)]))
        return unordered

    def _get_connected_wgraphs_w0(self, nvertices, nloops, nhairs):
        """Produces connected w-graphs with zero omega hairs. 
        Disregards the ordering of numbered hairs.
        Note that the graphs produced still have a (zero-valent) omega vertex,
        and a zero-or higher valent eps vertex.
        nloops is 2 minus the Euler characteristic of the returned graphs 
        """
        # print("wk0 ", nvertices, nloops, nhairs)
        # Algorithm: produce all hairy graphs, then fuse hairs
        max_eps_hairs = nloops+1
        for n_eps_hairs in range(max_eps_hairs+1):
            for G in self.get_hairy_graphs(nvertices, 1+nloops-n_eps_hairs, nhairs+n_eps_hairs, include_novertgraph=True):
                # make space for eps and omega vertex, and create both vertices
                # Warning: we assume here that we can modify the parameter... maybe make a copy first, but this costs space
                G.relabel({nvertices+j:nvertices+j+2 for j in range(nhairs+n_eps_hairs) })
                G.add_vertex(nvertices) # eps vertex
                G.add_vertex(nvertices+1) # omega vertex
                for S in itertools.combinations(range(nvertices+2, nvertices+nhairs+n_eps_hairs+2), n_eps_hairs):
                    GG = copy(G)
                    GG.merge_vertices([nvertices] + list(S))
                    if GG.size() < G.size():
                        continue
                    GG.relabel(range(G.order())) # make sure vertices are consecutively labeled
                    assert GG.order() == nvertices + nhairs + 2 , "Error: wrong graph size."
                    assert -GG.order() + 2 + GG.size() == nloops, "Error: wrong loop number"
                    yield GG

    def _get_connected_wgraphs_w0_memo(self, nvertices, nloops, nhairs):
        fName = "cache/get_connected_wgraphs_w0_%d_%d_%d.g6" % (nvertices, nloops, nhairs)
        part =  [ list(range(nvertices)) , [nvertices], [nvertices+1], list(range(nvertices+2, nvertices+nhairs+2))] 
        if os.path.isfile(fName):
            lst = StoreLoad.load_string_list(fName)
        else:
            lst = { G.canonical_label(partition=part).graph6_string() for G in self._get_connected_wgraphs_w0(nvertices, nloops, nhairs) }
            StoreLoad.store_string_list(lst, fName)

        return map(Graph, lst)

    def _get_connected_wgraphs_wk(self, nvertices, nloops, nhairs, nws):
        """Produces connected graphs with k omega-hairs. disregards the ordering of numbered hairs.
        """
        # print("wk ", nvertices, nloops, nhairs, nws)
        # Algorithm: produce all w graphs with nws hairs more, then fuse
        for G in self._get_connected_wgraphs_w0_memo(nvertices, nloops-nws, nhairs+nws):
            for S in itertools.combinations(range(nvertices+2, nvertices+nhairs+2+nws), nws):
                GG = copy(G)
                to_merge = [nvertices+1] + list(S)
                # print(to_merge, GG.vertices())
                GG.merge_vertices(to_merge)
                if GG.size() < G.size():
                    continue
                GG.relabel(range(G.order()))
                if not GG.has_edge(nvertices, nvertices+1):
                    yield GG
    
    def _get_connected_wgraphs_wk_memo(self, nvertices, nloops, nhairs, nws):
        fName = "cache/get_connected_wgraphs_wk_%d_%d_%d_%d.g6" % (nvertices, nloops, nhairs, nws)
        part =  [ list(range(nvertices)) , [nvertices], [nvertices+1], list(range(nvertices+2, nvertices+nhairs+2))] 
        if os.path.isfile(fName):
            lst = StoreLoad.load_string_list(fName)
        else:
            lst = { G.canonical_label(partition=part).graph6_string() for G in self._get_connected_wgraphs_wk(nvertices, nloops, nhairs, nws) }
            StoreLoad.store_string_list(lst, fName)

        return map(Graph, lst)


    def _get_pre_wgraphs_wk(self, nvertices, nloops, nhairs, nws):
        """Produces all wgraphs with at least one omega per connected component.
        Only implemented for nws =1,2.
        Disregards hair labels.
        """
        # print("pre wk ", nvertices, nloops, nhairs, nws)
        if nws == 1:
            # just forward the other routine
            for G in self._get_connected_wgraphs_wk_memo(nvertices, nloops, nhairs, 1):
                yield G
        elif nws == 2:
            # Our graph can have either one or two connected components
            # (A) one connected component.
            for G in self._get_connected_wgraphs_wk_memo(nvertices, nloops, nhairs, 2):
                yield G
            # (B) two connected components. 
            # We need to distribute hairs, loops and vertices over the components, in all possible ways
            # precompute graph lists

            for nvertices1 in range(int(nvertices/2)+1):
                for nloops1 in range(nloops+1):
                    for nhairs1 in range(nhairs+1):
                        # use itertools.product to avoid recomputing one list multiple times
                        for G1, G2 in itertools.product(self._get_connected_wgraphs_wk_memo(nvertices1, nloops1, nhairs1, 1),
                                            self._get_connected_wgraphs_wk_memo(nvertices-nvertices1, nloops-nloops1, nhairs-nhairs1, 1) ):
                            # take union of graphs, after relabeling
                            GG1 = G1.relabel( list(range(nvertices1))+list(range(nvertices, nvertices+2+nhairs1)), inplace=False)
                            GG2 = G2.relabel( list(range(nvertices1, nvertices+2))+list(range(nvertices+2+nhairs1, nvertices+2+nhairs)), inplace=False)
                            yield GG1.union(GG2)
        else:
            raise ValueError("Only nws=1 and 2 are implemented so far")

    def _get_pre_wgraphs_wk_memo(self, nvertices, nloops, nhairs, nws):
        fName = "cache/get_pre_wgraphs_wk_%d_%d_%d_%d.g6" % (nvertices, nloops, nhairs, nws)
        part =  [ list(range(nvertices)) , [nvertices], [nvertices+1], list(range(nvertices+2, nvertices+nhairs+2))] 
        if os.path.isfile(fName):
            lst = StoreLoad.load_string_list(fName)
        else:
            lst = { G.canonical_label(partition=part).graph6_string() for G in self._get_pre_wgraphs_wk(nvertices, nloops, nhairs, nws) }
            StoreLoad.store_string_list(lst, fName)

        return map(Graph, lst)

    def _attach_epsj(self,G, nvertices):
        """Takes the union of the w-graph G with an extra edge eps-hair
        nvertices must be the number of internal vertices, so that nvertices is also the label of the eps-vertex
        """
        GG = copy(G)
        n = GG.order()
        GG.add_vertex(n)
        GG.add_edge(nvertices, n)
        return GG
    
    def _get_all_wgraphs(self, nvertices, nloops, nhairs, nws):
        """Produces all (reduced) wgraphs.
        Disregards hair labels.
        """
        # Idea: just take the wgraphs produced above, and attach an extra hair and or an eps tadpole
        for G in self._get_pre_wgraphs_wk_memo(nvertices, nloops, nhairs, nws):
            yield G
        for G in self._get_pre_wgraphs_wk_memo(nvertices, nloops, nhairs-1, nws):
            yield self._attach_epsj(G, nvertices)
        # with tadpole
        for G in self._get_pre_wgraphs_wk_memo(nvertices, nloops-1, nhairs, nws):
            yield G
        for G in self._get_pre_wgraphs_wk_memo(nvertices, nloops-1, nhairs-1, nws):    
            yield self._attach_epsj(G, nvertices)

    def get_generating_graphs(self):
        # The routines above produce all wgraphs, we just have to permute the hair labels 
        
        # Produce all permutations of the hairs
        all_perm = [ list(range(0,self.n_vertices+2)) + list(p) for p in itertools.permutations(range(self.n_vertices+2, self.n_vertices+self.n_hairs+2)) ]

        return [G.relabel(p, inplace=False) for G in self._get_all_wgraphs(self.n_vertices, self.n_loops, self.n_hairs, self.n_ws) 
                                            for p in all_perm ]

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
            if n_l == 1: #hair
                continue
            elif n_l == 2: #edge
                G.add_edge(neighbors)
                G.delete_vertex(v)
            else:
                raise ValueError('%s: Vertices of second colour should have 1 or 2 neighbours' % str(self))
        G.relabel(range(0,G.order()))
        return G


class WRHairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of hairy graph vector spaces with specified number of omega hairs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_range (range): Range for the number of hairs.
        - w_range (range): number of omega hairs
        - sub_type (str): Sub type of graphs.
    """
    def __init__(self, v_range, l_range, h_range, w_range):
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
        self.w_range = w_range
        self.sub_type = ""

        vs_list = [WRHairyGraphVS(v, l, h, w) for
                   (v, l, h, w) in itertools.product(self.v_range, self.l_range, self.h_range, self.w_range)]
        super(WRHairyGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return 'wrhairy graphs'

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('hairs', self.h_range), ('ws', self.w_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s" % graph_type
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)


# ------- Operators --------
class ContractEdgesGO(GraphOperator.GraphOperator):
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
               and domain.n_hairs == target.n_hairs and domain.n_ws == target.n_ws

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, n_ws):
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
        domain = WRHairyGraphVS(n_vertices, n_loops, n_hairs, n_ws)
        target = WRHairyGraphVS(n_vertices - 1, n_loops, n_hairs, n_ws)
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
            # only edges not connected to a numbered hair-vertex can be contracted
            if u >= self.domain.n_vertices+2 or v >= self.domain.n_vertices+2:
                continue

            # ensure u<v (this should be always true anyway actually)
            if u>v:
                u,v = v,u
            
            sgn = 1 if i % 2 == 0 else -1
            previous_size = G.size()
            previous_has_tadpole = (previous_size - self.domain.n_vertices - self.domain.n_hairs < self.domain.n_loops) 
            sgn *= -1 if previous_has_tadpole else 1
            G1 = copy(G)
            # label all edges to determine sign later
            Shared.enumerate_edges(G1)

            # we always delete the lower index vertex. This ensures that the extra vertices are never deleted
            if v <= self.domain.n_vertices:
                G1.merge_vertices([v, u])
                if (previous_size - G1.size()) != 1:
                    continue
                G1.relabel(range(0,self.domain.n_vertices+1+self.domain.n_hairs), inplace=True)
                # find edge permutation sign 
                sgn *= Shared.shifted_edge_perm_sign2(G1)
                image.append((G1, sgn))
            elif v==self.domain.n_vertices+1:
                # the second vertex is now omega, so we need to merge the vertex with the eps vertex
                # after reonnecting one of the edges to omega
                # we assume that u != eps, because this is forbidden in our graphs
                G1.delete_edge(u,v)
                # special care must be taken since a tadpole could be created at eps
                # and this is true iff there is an edge u-eps
                eps = self.domain.n_vertices
                new_has_tadpole = G1.has_edge(u,eps)
                # double tadpole => zero
                if new_has_tadpole and previous_has_tadpole:
                    continue
                if new_has_tadpole:
                    # remove the edge and compute the appropriate sign
                    k = G1.edge_label(u,eps)
                    G1.delete_edge(u,eps)
                    sgn *= 1 if ( (k % 2==0) == (k<i) ) else -1
                # loop over other neighbors w to be connected to omega
                for w in G1.neighbors(u):
                    G2 = copy(G1)
                    sgn2 = sgn
                    # reconnect the w-v-edge to omega (i.e., to v)
                    old_label = G2.edge_label(u,w)
                    G2.delete_edge(u,w)
                    G2.add_edge(w,v, old_label)
                    
                    # now merge u and eps
                    G2.merge_vertices([eps, u])
                    # in case we have too few edges some double edges have been created => zero
                    if (previous_size - G2.size()) != (2 if new_has_tadpole else 1):
                        continue
                    G2.relabel(range(0,self.domain.n_vertices+1+self.domain.n_hairs), inplace=True)
                    # find edge permutation sign 
                    sgn2 *= Shared.shifted_edge_perm_sign2(G2)
                    image.append((G2, sgn2))

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

# ------- Graph Complex --------
class WRHairyGC(GraphComplex.GraphComplex):
    """Graph complex for hairy graphs.

    Attributes:
        - v_range (range): Range for the number of vertices.
        - l_range (range): Range for the number of loops.
        - h_range (range): Range for the number of hairs.
        - w_range (range): Range for the number of omega hairs
        - sub_type (str): Sub type of graphs.
    """
    def __init__(self, v_range, l_range, h_range, w_range, differentials):
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
        self.w_range = w_range
        self.sub_type = ""

        sum_vector_space = WRHairyGraphSumVS(self.v_range, self.l_range, self.h_range, self.w_range)
        differential_list = []
        if not set(differentials).issubset(['contract']):
            raise ValueError("Differentials for hairy graph complex: 'contract'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)
        super(WRHairyGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))

    def print_dim_and_eulerchar(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    ds = [WRHairyGraphVS(v,l,h,w).get_dimension() for v in self.v_range ]
                    eul = sum( [(1 if j%2==0 else -1) * d for j,d in enumerate(ds)] )
                    print("Dimensions (w,h,l) ",w,h,l, ":", ds, "Euler", eul)
    def print_cohomology_dim(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    cohomdict = {}
                    for v in self.v_range:
                        D1 = ContractEdgesGO.generate_operator(v,l,h,w)
                        D2 = ContractEdgesGO.generate_operator(v+1,l,h,w)
                        try:
                            d = WRHairyGraphVS(v,l,h,w).get_dimension()
                            r1 = D1.get_matrix_rank()
                            r2 = D2.get_matrix_rank()
                            cohomdict[v] = d-r1-r2
                        except:
                            pass
                    
                    print("Cohomology Dimensions (w,h,l) ",w,h,l, ":", cohomdict)