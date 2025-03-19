"""Graph complexes based on simple graphs with numbered hairs and hairs of two (omega- and epsilon-)decorations 
as in the graph complex computing weight 11 cohomology.

The omega decorations are considered odd !!!!

Implemented Differentials: Contract edges, make eps vertex to omega vertex.


Graphs are realized as simple graphs with 1+w extra vertices for epsilon and omega, (index self.n_vertices and self.n_vertices+1,2,...,w).
WARNING: If there is an eps-eps tadpole the corresponding loop is not part of the graph--one can determine the presence of the tadpole
from the overall one too small loop number.
There must not be tadpoles at internal vertices
TODO: Take care that this does not produce problems
"""


# __all__ = ['WOHairyGraphVS', 'WOHairyGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
#            'RestrictedContractEdgesGO', 'RestrictedContractEdgesD',
#            'SymmProjector', 'WOHairyGC']

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
import inspect
import CHairyGraphComplex
from sage.combinat.shuffle import ShuffleProduct

graph_type = "wohairy2"


# unshuffles_dict = {}

# def unshuffles_old(p,q):
#     """
#     Returns all p,q-unshuffles. Thes are all permutations [i1,...ipq] of 0,..,p+q-1 such that i1<..<ip and ip+1<..<ipq
#     """
#     if (p,q) in unshuffles_dict:
#         return unshuffles_dict[(p,q)]
#     else:
#         ret = [[j for j in pp] for pp in itertools.permutations(range(p+q)) 
#             if all(pp[i] < pp[i+1] for i in range(p-1))
#             and all(pp[i] < pp[i+1] for i in range(p,p+q-1))]
#         unshuffles_dict[(p,q)] = ret
#         return ret
    
def unshuffles(p,q):
    for pp in itertools.combinations(range(p+q), p):
        yield list(pp) + [j for j in range(p+q) if j not in pp]

def dump_args(func):
    """
    Decorator to print function call details.

    This includes parameters names and effective values.
    """

    def wrapper(*args, **kwargs):
        func_args = inspect.signature(func).bind(*args, **kwargs).arguments
        func_args_str = ", ".join(
            map("{0[0]} = {0[1]!r}".format, func_args.items()))
        print(f"{func.__module__}.{func.__qualname__} ( {func_args_str} )")
        return func(*args, **kwargs)

    return wrapper

# ------- Graph Vector Space --------




class WOHairyGraphPreVS(SymmetricGraphComplex.SymmetricGraphVectorSpace):
    """WOHairy graph vector space.

    Sub vector space with specified number of (internal) vertices, loops, numbered hairs, and omega decorations (ws)
    and at least trivalent vertices. One hair is composed of a hair vertex and an edge connecting it to a vertex.

    Attributes:
        - n_vertices (int): Number of internal vertices.
        - n_loops (int): Number of loops.
        - n_hairs (int): Number of hairs.
        - n_edges (int): Number of edges, not counting edges to numbered hairs.
        - sub_type (str): Sub type of graphs. This is currently not used, but might be needed for an extension later.
        - ogvs (OrdinaryGraphComplex.OrdinaryGVS): Ordinary graph vector space without hairs.

    """

    def __init__(self, n_vertices, n_loops, n_hairs, n_ws, n_comp):
        """Initialize the hairy graph vector space.

        :param n_vertices: Number of internal vertices.
        :type n_vertices: int
        :param n_loops: the genus of the graph. Warning: this is counted without the special vertex.... e.g., in weight 11 the actual genus is one more
        :type n_loops: int
        :param n_hairs: Number of hairs. They are distinguishable, numbered 1,...,n
        :type n_hairs: int
        :param n_ws: Number of omega decorated hairs.
        :type n_ws: int
        :param n_comp: Number of connected components.
        :type n_comp: int
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs = n_hairs
        self.n_ws = n_ws
        self.n_comp = n_comp

        # we count only the internal edges and omega and eps edges, but not the hair edges
        self.n_edges = self.n_loops + self.n_vertices
        self.sub_type = ""
        super(WOHairyGraphPreVS, self).__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(
            self.n_vertices + self.n_hairs+1+self.n_ws, self.n_loops, False)

    def get_type(self):
        return graph_type

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs and self.n_ws == other.n_ws and self.n_comp == other.n_comp

    def __hash__(self):
        return hash("wograpre%d_%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wograpre%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "wograpre%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('hairs', self.n_hairs), ('ws', self.n_ws), ('comp', self.n_comp)])

    def get_partition(self):
        # All internal vertices are in color 1, the single eps vertex in color 2, the w vertices in color 3
        # and the hair vertices are in colors 4,...,n+3.
        ret = [list(range(0, self.n_vertices))] + [[self.n_vertices]] + [list(range(self.n_vertices+1, self.n_vertices +
                                                                                    1+self.n_ws))] + [[j] for j in range(self.n_vertices+1+self.n_ws, self.n_vertices + self.n_hairs+self.n_ws+1)]
        # take out empty lists
        # return [l for l in ret if l != []]
        return ret

    def plot_graph(self, G):
        GG = Graph(G, loops=True)
        # add proper tadpole if needed
        if GG.size() < self.n_edges + self.n_hairs:
            GG.add_edge(self.n_vertices, self.n_vertices)

        return GG.plot(partition=self.get_partition(), vertex_labels=True)

    def is_valid(self):
        # At least trivalent internal vertices.
        l = (3 * self.n_vertices + self.n_ws <=
             2 * self.n_edges + self.n_hairs)
        # Nonnegative number of vertices, non negative number of loops, non-negative number of hairs,
        l = l and self.n_vertices >= 0 and self.n_loops >= 0 and self.n_hairs >= 0 and self.n_edges >= 0 and self.n_comp >=1
        # At most a full graph.
        # l = l and self.n_edges <= (
        #     self.n_vertices+2) * (self.n_vertices + 1) / 2
        
        l = l and (self.n_ws <= self.n_loops + self.n_comp)
        # excess must be non-negative
        l = l and (2*self.n_ws <= 3*self.n_loops+2*self.n_hairs)
        # each connected component must contribute to the genus, or otherwise have a hair
        l = l and (self.n_comp <= self.n_hairs + self.n_loops)
        return l

    def get_work_estimate(self):
        # TODO
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return GCDimensions.get_wrhairy_dim_estimate(self.n_vertices, self.n_loops, self.n_hairs, self.n_ws)
        # return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) * (self.n_vertices ** self.n_hairs) / factorial(self.n_vertices)

    def get_chairy_prerequisites(self):
        if self.n_comp == 1:
            return [CHairyGraphComplex.CHairyGraphVS(self.n_vertices, l, self.n_hairs+self.n_loops+1-l, False) 
                     for l in range(self.n_loops - self.n_ws +2) ]
        else:
            ret = [CV for V in self.get_self_prerequisites() for CV in V.get_chairy_prerequisites() ]
            return list(set(ret))
    def get_self_prerequisites(self):
        if self.n_comp == 1:
            return []
        else:
            ret = []
            for l in range(self.n_loops+1):
                for v in range(self.n_vertices+1):
                    for h in range(self.n_hairs+1):
                        for w in range(self.n_ws+1):
                            V1 = WOHairyGraphPreVS(v, l, h, w, self.n_comp-1)
                            V2 = WOHairyGraphPreVS(self.n_vertices-v, self.n_loops-l, self.n_hairs-h, self.n_ws-w, 1)
                            if V1.is_valid() and V2.is_valid():
                                ret += [V1, V2]
                                ret += V1.get_self_prerequisites()
            return list(set(ret))

    def get_unshuffles(self, p,q, offset=0):
        """
        Returns all p,q-unshuffles. Thes are all permutations [i1,...ipq] of 0,..,p+q-1 such that i1<..<ip and ip+1<..<ipq
        """
        return [[j+offset for j in pp] for pp in unshuffles(p,q)]

    def get_generating_graphs(self):
        """
        Produces a list of generating graphs
        """
        nvertices = self.n_vertices
        nloops = self.n_loops
        nhairs = self.n_hairs
        nws = self.n_ws

        if self.n_comp == 1:
            # case with no vertices has to be treated separately, since not contained by default in CHairyGraphVS
            if nvertices == 0:
                if nws==0 and nhairs == 0 and nloops ==1:
                    #yield Graph(1) # eps-eps edge
                    # for now we forbid the eps-eps edge. It has to be added later
                    pass
                elif (nws==1 and nhairs == 0 and nloops ==1): # eps-omega 
                    G = Graph(2)
                    G.add_edge(0,1)
                    yield G
                elif (nws==0 and nhairs == 1 and nloops ==0): # eps-j
                    G = Graph(2)
                    G.add_edge(0,1)
                    yield G
                elif (nws==1 and nhairs == 1 and nloops ==0): # omega-j
                    G = Graph(3)
                    G.add_edge(1,2)
                    yield G
                return

            for l in range(self.n_loops - self.n_ws +2):
                neps = nloops - nws +1 - l
                V = CHairyGraphComplex.CHairyGraphVS(nvertices, l, nhairs+nws+neps, False) 
                if not V.is_valid():
                    continue
                # fuse eps vertices. Add one first to cover case eps=0
                for G in V.get_basis():
                    # G = copy(tG)
                    # print(G.graph6_string())
                    oldsize = G.size()
                    G.add_vertex()
                    G.merge_vertices(list(range(nvertices,nvertices+neps ))+[nvertices+nws+nhairs+neps])
                    # need to relabel so that vertices are consecutively labeled after merge
                    G.relabel()
                    # print("merged", G.graph6_string())
                    # print(G.size(), self.n_edges)
                    #print(self.n_edges)
                    # if G has too few edges we had a double eps edge-> graph is zero
                    if G.size() == oldsize:
                        yield G
        else: # n_comp >1
            for l in range(self.n_loops+1):
                for v in range(self.n_vertices+1):
                    for h in range(self.n_hairs+1):
                        for w in range(self.n_ws+1):
                            V1 = WOHairyGraphPreVS(v, l, h, w, self.n_comp-1)
                            V2 = WOHairyGraphPreVS(self.n_vertices-v, self.n_loops-l, self.n_hairs-h, self.n_ws-w, 1)
                            if V1.is_valid() and V2.is_valid():
                                for G1 in V1.get_basis():
                                    for G2 in V2.get_basis():
                                        GG = G1.disjoint_union(G2)
                                        GG.relabel()
                                        # merge the eps vertices 
                                        GG.merge_vertices([v, nvertices+1+w+h])
                                        # reorder so that vertices are aligned
                                        GG.relabel(list(range(v)) # old internal vertices 1
                                                    + [nvertices] # new eps vertex 
                                                    + list(range(nvertices+1, nvertices+w+1)) # old ws 1
                                                    + list(range(nvertices+nws+1, nvertices+nws+1+h)) # old hairs 1
                                                    + list(range(v,nvertices)) # old internal vertices 2
                                                    + list(range(nvertices+1+w, nvertices+nws+1)) # old ws 2
                                                    + list(range(nvertices+1+nws+h,nvertices+1+nws+nhairs)) # old hairs 2
                                                    )
                                        # sum over all shuffle permutations, apply to the hairs
                                        # list of all h, nhairs-h - shuffles:
                                        vlist = list(range(nvertices+1+nws))
                                        all_shuffle = [vlist + p for p in self.get_unshuffles(h,nhairs-h, nvertices+1+nws)]
                                        for sh in all_shuffle:
                                            GGG = GG.relabel(sh, inplace=False)
                                            yield GGG


    # def get_generating_graphs(self):
    #     # The routines above produce all wgraphs, we just have to permute the hair labels

    #     # Produce all permutations of the hairs
    #     all_perm = [list(range(0, self.n_vertices+2)) + list(p)
    #                 for p in itertools.permutations(range(self.n_vertices+2, self.n_vertices+self.n_hairs+2))]

    #     return (G.relabel(p, inplace=False) for G in self._get_all_wgraphs(self.n_vertices, self.n_loops, self.n_hairs, self.n_ws)
    #             for p in all_perm)

    def perm_sign(self, G, p):
        # The sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the omega-vertices.
        sgn = self.ogvs.perm_sign(G, p)

        # Compute the extra contribution from omegas.
        # if self.even_hairs == self.even_edges:
        if self.n_ws > 0:
            wperm = p[self.n_vertices+1:]
            if len(wperm) != 0:
                sgn *= Shared.Perm.shifted(wperm).signature()
        return sgn

    def get_n(self):
        return self.n_hairs

    def vertex_permutation_from_permutation(self, p):
        return list(range(0, self.n_vertices+2)) + [j+self.n_vertices+1 for j in p]

    def get_isotypical_projector(self, rep_index):
        return SymmProjector(self, rep_index)




class WOHairyGraphVS(SymmetricGraphComplex.SymmetricGraphVectorSpace):
    """WOHairy graph vector space.

    Sub vector space with specified number of (internal) vertices, loops, numbered hairs, and omega decorations (ws)
    and at least trivalent vertices. One hair is composed of a hair vertex and an edge connecting it to a vertex.

    Attributes:
        - n_vertices (int): Number of internal vertices.
        - n_loops (int): Number of loops.
        - n_hairs (int): Number of hairs.
        - n_edges (int): Number of edges, not counting edges to numbered hairs.
        - sub_type (str): Sub type of graphs. This is currently not used, but might be needed for an extension later.
        - ogvs (OrdinaryGraphComplex.OrdinaryGVS): Ordinary graph vector space without hairs.

    """

    def __init__(self, n_vertices, n_loops, n_hairs, n_ws):
        """Initialize the hairy graph vector space.

        :param n_vertices: Number of internal vertices.
        :type n_vertices: int
        :param n_loops: the genus of the graph. Warning: this is counted without the special vertex.... e.g., in weight 11 the actual genus is one more
        :type n_loops: int
        :param n_hairs: Number of hairs. They are distinguishable, numbered 1,...,n
        :type n_hairs: int
        :param n_ws: Number of omega decorated hairs.
        :type even_edges: int
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.n_hairs = n_hairs
        self.n_ws = n_ws
        # we count only the internal edges and omega and eps edges, but not the hair edges
        self.n_edges = self.n_loops + self.n_vertices
        self.sub_type = ""
        super(WOHairyGraphVS, self).__init__()
        self.ogvs = OrdinaryGraphComplex.OrdinaryGVS(
            self.n_vertices + self.n_hairs+1+self.n_ws, self.n_loops, False)

    def get_type(self):
        return graph_type

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops and self.n_hairs == other.n_hairs and self.n_ws == other.n_ws

    def __hash__(self):
        return hash("wogra%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wogra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_basis_file_path(self):
        s = "wogra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('hairs', self.n_hairs), ('ws', self.n_ws)])

    def get_partition(self):
        # All internal vertices are in color 1, the single eps vertex in color 2, the w vertex in color 3
        # and the hair vertices are in colors 4,...,n+3.
        ret = [list(range(0, self.n_vertices))] + [[self.n_vertices]] + [list(range(self.n_vertices+1, self.n_vertices +
                                                                                    1+self.n_ws))] + [[j] for j in range(self.n_vertices+1+self.n_ws, self.n_vertices + self.n_hairs+self.n_ws+1)]
        # take out empty lists
        # return [l for l in ret if l != []]
        return ret

    def plot_graph(self, G):
        GG = Graph(G, loops=True)
        # add proper tadpole if needed
        if GG.size() < self.n_edges + self.n_hairs:
            GG.add_edge(self.n_vertices, self.n_vertices)

        return GG.plot(partition=self.get_partition(), vertex_labels=True)

    def is_valid(self):
        # At least trivalent internal vertices.
        l = (3 * self.n_vertices + self.n_ws <=
             2 * self.n_edges + self.n_hairs)
        # Nonnegative number of vertices, non negative number of loops, non-negative number of hairs,
        l = l and self.n_vertices >= 0 and self.n_loops >= 0 and self.n_hairs >= 0 and self.n_edges >= 0
        # At most a full graph.
        # l = l and self.n_edges <= (
        #     self.n_vertices+2) * (self.n_vertices + 1) / 2
        return l

    def get_work_estimate(self):
        # TODO
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return GCDimensions.get_wrhairy_dim_estimate(self.n_vertices, self.n_loops, self.n_hairs, self.n_ws)
        # return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) * (self.n_vertices ** self.n_hairs) / factorial(self.n_vertices)

    def get_prerequisites(self):
        # need to add all graphs with same parameters as myself, or with one edge lessbecause of eps-eps-edge
        maxc = self.n_ws+self.n_loops+self.n_hairs # crude estimate
        ret = [WOHairyGraphPreVS(self.n_vertices,self.n_loops, self.n_hairs, self.n_ws, c) for c in range(maxc+1)] + \
            [WOHairyGraphPreVS(self.n_vertices,self.n_loops-1, self.n_hairs, self.n_ws, c) for c in range(maxc+1)]
        ret = [V for V in ret if V.is_valid()]
        return ret


    def get_generating_graphs(self):
        """
        Produces a list of generating graphs
        """
        for V in self.get_prerequisites():
            for G in V.get_basis():
                yield G

    # def get_generating_graphs(self):
    #     # The routines above produce all wgraphs, we just have to permute the hair labels

    #     # Produce all permutations of the hairs
    #     all_perm = [list(range(0, self.n_vertices+2)) + list(p)
    #                 for p in itertools.permutations(range(self.n_vertices+2, self.n_vertices+self.n_hairs+2))]

    #     return (G.relabel(p, inplace=False) for G in self._get_all_wgraphs(self.n_vertices, self.n_loops, self.n_hairs, self.n_ws)
    #             for p in all_perm)

    def perm_sign(self, G, p):
        # The sign is the same as the corresponding sign in the
        # ordinary graph complex, apart from an extra contribution from the omega-vertices.
        sgn = self.ogvs.perm_sign(G, p)

        # Compute the extra contribution from omegas.
        # if self.even_hairs == self.even_edges:
        if self.n_ws > 0:
            wperm = p[self.n_vertices+1:]
            if len(wperm) != 0:
                sgn *= Shared.Perm.shifted(wperm).signature()
        return sgn

    def get_n(self):
        return self.n_hairs

    def vertex_permutation_from_permutation(self, p):
        return list(range(0, self.n_vertices+2)) + [j+self.n_vertices+1 for j in p]

    def get_isotypical_projector(self, rep_index):
        return SymmProjector(self, rep_index)


class WOHairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
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

        vs_list = [WOHairyGraphVS(v, l, h, w) for
                   (v, l, h, w) in itertools.product(self.v_range, self.l_range, self.h_range, self.w_range)]
        super(WOHairyGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return 'wohairy graphs'

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range), ('hairs', self.h_range), ('ws', self.w_range)])

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
            and domain.n_hairs == target.n_hairs and domain.n_ws == target.n_ws

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, n_ws):
        """Return a contract edges graph operator.

        :param n_vertices: Number of vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param n_hairs: Number of (numbered) hairs.
        :type n_hairs: int
        :param n_ws: Number of omega hairs.
        :type n_ws: int
        :return: Contract edges graph operator based on the specified domain vector space.
        :rtype: ContractEdgesGO
        """
        domain = WOHairyGraphVS(n_vertices, n_loops, n_hairs, n_ws)
        target = WOHairyGraphVS(n_vertices - 1, n_loops, n_hairs, n_ws)
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
        # Returns as work estimate: domain.n_edges * domain_dim * log(target dimension, 2)
        if not self.is_valid():
            return 0
        try:
            (domain_dim, dimtarget_dim) = (
                self.domain.get_dimension(), self.target.get_dimension())
        except StoreLoad.FileNotFoundError:
            return 0
        if domain_dim == 0 or dimtarget_dim == 0:
            return 0
        return self.domain.n_edges * domain_dim * math.log(dimtarget_dim, 2)

    def get_type(self):
        return 'contract edges'

    def operate_on(self, G):
        # print("operate on:", G.graph6_string(),
            #   self.domain.get_ordered_param_dict())
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image = []
        for (i, e) in enumerate(G.edges(labels=False,sort=True)):
            (u, v) = e

            # ensure u<v (this should be always true anyway actually)
            if u > v:
                u, v = v, u

            # only edges connected to at least one internal vertex, and not connected to a numbered hair-vertex can be contracted
            if u >= self.domain.n_vertices or v >= self.domain.n_vertices+self.domain.n_ws+1:
                continue

            sgn = 1 if i % 2 == 0 else -1
            # print("sgn0",sgn)
            previous_size = G.size()
            previous_has_tadpole = (
                previous_size - self.domain.n_vertices - self.domain.n_hairs < self.domain.n_loops)
            sgn *= -1 if previous_has_tadpole else 1
            # print("sgn1",sgn)
            G1 = copy(G)
            # label all edges to determine sign later
            Shared.enumerate_edges(G1)

            # we always delete the lower index vertex. This ensures that the extra vertices are never deleted
            if v <= self.domain.n_vertices:
                G1.merge_vertices([v, u])
                if (previous_size - G1.size()) != 1:
                    continue
                G1.relabel(range(0, self.domain.n_vertices+self.domain.n_ws +
                           self.domain.n_hairs), inplace=True)
                # find edge permutation sign
                sgn *= Shared.shifted_edge_perm_sign2(G1)
                # print("sgn3_",sgn)
                image.append((G1, sgn))
                # image.append((Graph(G1.graph6_string()), sgn))
                # print("hmm0:", G.graph6_string(), G1.graph6_string())
            elif u < self.domain.n_vertices and v >= self.domain.n_vertices+1:
                # the second vertex is now an omega-vertex, so we need to merge the vertex with the eps vertex
                # after reconnecting one of the edges to omega
                # we assume that u != eps, because eps-omega-edges cannot be contracted
                G1.delete_edge(u, v)
                # special care must be taken since a tadpole could be created at eps
                # and this is true iff there is an edge u-eps
                eps = self.domain.n_vertices
                # new_has_tadpole = G1.has_edge(u, eps)
                # # double tadpole => zero
                # if new_has_tadpole and previous_has_tadpole:
                #     continue
                # if new_has_tadpole:
                #     # remove the edge and compute the appropriate sign
                #     k = G1.edge_label(u, eps)
                #     G1.delete_edge(u, eps)
                #     sgn *= 1 if ((k % 2 == 0) == (k < i)) else -1

                # loop over neighbors w to be connected to omega
                for w in G1.neighbors(u):
                    G2 = copy(G1)
                    sgn2 = sgn
                    # reconnect the w-v-edge to omega (i.e., to v)
                    old_label = G2.edge_label(u, w)
                    G2.delete_edge(u, w)
                    G2.add_edge(w, v, old_label)

                    # we want to merge u and eps... however, this might create a tadpole
                    new_has_tadpole = G2.has_edge(u, eps)
                    if new_has_tadpole and previous_has_tadpole:
                        continue
                    if new_has_tadpole:
                        # remove the edge and compute the appropriate sign
                        k = G2.edge_label(u, eps)
                        G2.delete_edge(u, eps)
                        sgn2 *= 1 if ((k % 2 == 0) == (k < i)) else -1

                    # now merge u and eps
                    G2.merge_vertices([eps, u])
                    # in case we have too few edges some double edges have been created => zero
                    if (previous_size - G2.size()) != (2 if new_has_tadpole else 1):
                        continue
                    G2.relabel(range(0, self.domain.n_vertices +
                               self.domain.n_hairs+self.domain.n_ws), inplace=True)
                    # find edge permutation sign
                    sgn2 *= Shared.shifted_edge_perm_sign2(G2)
                    # sanity checks
                    if G2.order() != self.target.n_vertices+self.target.n_hairs+self.target.n_ws+1:
                        print("Error contract:", G.graph6_string(),
                              G2.graph6_string())
                    # else:
                    #     print("hmm:", G.graph6_string(), G2.graph6_string())
                    image.append((G2, sgn2))

        return image

    def restrict_to_isotypical_component(self, rep_index):
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
        s = "cohomology_dim_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_D_%s" % (graph_type)
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


class EpsToOmegaGO(SymmetricGraphComplex.SymmetricGraphOperator):
    """Operator that makes one eps into an omega hair.

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
        super(EpsToOmegaGO, self).__init__(domain, target)

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
        return domain.n_vertices == target.n_vertices and domain.n_loops == target.n_loops \
            and domain.n_hairs == target.n_hairs and domain.n_ws+1 == target.n_ws

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, n_hairs, n_ws):
        """Return a contract edges graph operator.

        :param n_vertices: Number of vertices of the domain.
        :type n_vertices: int
        :param n_loops: Number of loops of the domain.
        :type n_loops: int
        :param n_hairs: Number of (numbered) hairs.
        :type n_hairs: int
        :param n_ws: Number of omega hairs.
        :type n_ws: int
        :return: Contract edges graph operator based on the specified domain vector space.
        :rtype: ContractEdgesGO
        """
        domain = WOHairyGraphVS(n_vertices, n_loops, n_hairs, n_ws)
        target = WOHairyGraphVS(n_vertices, n_loops, n_hairs, n_ws+1)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "epstowD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "epstowD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    # def get_ref_matrix_file_path(self):
    #     s = "epstowD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    # def get_ref_rank_file_path(self):
    #     s = "epstowD%d_%d_%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
    #     return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * domain_dim * log(target dimension, 2)
        if not self.is_valid():
            return 0
        try:
            (domain_dim, dimtarget_dim) = (
                self.domain.get_dimension(), self.target.get_dimension())
        except StoreLoad.FileNotFoundError:
            return 0
        if domain_dim == 0 or dimtarget_dim == 0:
            return 0
        return self.domain.n_edges * domain_dim * math.log(dimtarget_dim, 2)

    def get_type(self):
        return 'contract edges'

    def operate_on(self, G):
        G1 = copy(G)
        sgn = (-1)**G.size()
        image = []

        previous_size = G.size()
        previous_has_tadpole = (
            previous_size - self.domain.n_vertices - self.domain.n_hairs < self.domain.n_loops)
        # if we have a tadpole, there is secretely one more edge that is not included in the graph
        if previous_has_tadpole:
            sgn *= -1

        # label all edges to determine sign later
        Shared.enumerate_edges(G1)

        # add one new omega vertex in position n_vertices +1
        G1.relabel(list(range(self.domain.n_vertices+1)) + list(range(self.domain.n_vertices +
                   2, self.domain.n_vertices+self.domain.n_ws+self.domain.n_hairs+2)))
        G1.add_vertex(self.domain.n_vertices+1)

        # reconnect one eps edge to the new vertex
        eps = self.domain.n_vertices
        new_w = self.domain.n_vertices+1
        for v in G1.neighbors(eps):
            sgn2 = sgn
            G2 = copy(G1)
            old_label = G2.edge_label(v, eps)
            G2.delete_edge(v, eps)
            G2.add_edge(v, new_w, old_label)
            sgn2 *= Shared.shifted_edge_perm_sign2(G2)
            image.append((G2, sgn2))

        # In case the original graph has a tadpole at eps, we also have to reconnect the tadpole edge (twice)
        if previous_has_tadpole:
            sgn2 = sgn
            G2 = copy(G1)
            # the tadpole edges label is the first in the ordering
            G2.add_edge(eps, new_w, -1)
            sgn2 *= Shared.shifted_edge_perm_sign2(G2)
            # factor 2 because tadpole has 2 eps vertices
            image.append((G2, 2*sgn2))

        return image

    def restrict_to_isotypical_component(self, rep_index):
        return RestrictedEpsToOmegaGO(self, rep_index)


class RestrictedEpsToOmegaGO(SymmetricGraphComplex.SymmetricRestrictedOperatorMatrix):
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
        return EpsToOmegaGO.is_match(domain.vs, target.vs) and domain.rep_index == target.rep_index


class EpsToOmegaD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: HairyGraphSumVS
        """
        super(EpsToOmegaD, self).__init__(sum_vector_space,
                                          EpsToOmegaGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'EpsToOmega'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class RestrictedEpsToOmegaD(SymmetricGraphComplex.SymmetricDifferential):

    def get_type(self):
        return 'isotypical epstoomega'

    def get_cohomology_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.diff.sum_vector_space.sub_type
        s = "info_epstoomega_D_iso_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class SymmProjector(SymmetricGraphComplex.SymmetricProjectionOperator):
    """This class encodes the projector to an isotypical component of the symmetric group action
        by permuting numbered hairs.
        Warning: The matrix stores not the projector, but projector * n_hairs! / rep_dimension??, to have integral matrices.

    Attributes:
        - sub_type(str): Graphs sub type of the domain.
    """

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
        return Shared.OrderedDict([('vertices', do.n_vertices), ('loops', do.n_loops), ('hairs', do.n_hairs), ('ws', do.n_ws), ('rep_index', self.rep_index)])

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


# ------- Graph Complex --------
class WOHairyGC(GraphComplex.GraphComplex):
    """Graph complex for hairy graphs.

    Attributes:
        - v_range(range): Range for the number of vertices.
        - l_range(range): Range for the number of loops.
        - h_range(range): Range for the number of hairs.
        - w_range(range): Range for the number of omega hairs
        - sub_type(str): Sub type of graphs.
    """

    def __init__(self, v_range, l_range, h_range, w_range, differentials):
        """Initialize the graph complex.

        : param v_range: Range for the number of vertices.
        : type v_range: range
        : param l_range: Range for the number of loops.
        : type l_range: range
        : param h_range: Range for the number of numbered hairs.
        : type  h_range: range
        : param w_range: Range for number of omega hairs.
        : type w_range: range
        : param differentials: List of differentials. Options: 'contract', 'et1h'.
        : type differentials: list(str)
        """
        self.v_range = v_range
        self.l_range = l_range
        self.h_range = h_range
        self.w_range = w_range
        self.sub_type = ""

        sum_vector_space = WOHairyGraphSumVS(
            self.v_range, self.l_range, self.h_range, self.w_range)
        differential_list = []
        if not set(differentials).issubset(['contract', 'contract_iso', 'epstoomega', 'epstoomega_iso']):
            raise ValueError(
                "Differentials for hairy graph complex: 'contract'")
        contract_edges_dif = ContractEdgesD(sum_vector_space)
        epstoomega_dif = EpsToOmegaD(sum_vector_space)
        if 'contract' in differentials:
            differential_list.append(contract_edges_dif)
        if 'contract_iso' in differentials:
            contract_iso_edges_dif = RestrictedContractEdgesD(
                contract_edges_dif)
            differential_list.append(contract_iso_edges_dif)
            print("Attention: contract_iso operates on nonzero cohomology entries only, so they need to be computed before!")
        if 'epstoomega' in differentials:
            differential_list.append(epstoomega_dif)
        super(WOHairyGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))

    def print_dim_and_eulerchar(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    ds = [WOHairyGraphVS(v, l, h, w).get_dimension()
                          for v in self.v_range]
                    eul = sum([(1 if j % 2 == 0 else -1) *
                              d for j, d in enumerate(ds)])
                    print("Dimensions (w,h,l) ", w,
                          h, l, ":", ds, "Euler", eul)

    def print_cohomology_dim(self):
        for w in self.w_range:
            for h in self.h_range:
                for l in self.l_range:
                    cohomdict = {}
                    for v in self.v_range:
                        D1 = ContractEdgesGO.generate_operator(v, l, h, w)
                        D2 = ContractEdgesGO.generate_operator(v+1, l, h, w)
                        try:
                            d = WOHairyGraphVS(v, l, h, w).get_dimension()
                            r1 = D1.get_matrix_rank()
                            r2 = D2.get_matrix_rank()
                            cohomdict[v] = d-r1-r2
                        except:
                            pass

                    print("Cohomology Dimensions (w,h,l) ",
                          w, h, l, ":", cohomdict)


# --------- Bicomplex ------------------

# class WOHairyDegSlice(SymmetricGraphComplex.SymmetricDegSlice):
#     """Degree slice of forested graphs

#     Total degree = marked edges

#     Attributes:
#         - even_edges (bool): True for even edges, False for odd edges.

#     """

#     def is_complete(self):
#         for vs in self.vs_list:
#             if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
#                 return False
#         return True

#     def __init__(self, n_loops, n_marked_edges, n_hairs, even_edges):
#         """Initialize the degree slice.

#         :param deg: Total degree of the degree slice.
#         :type deg: int
#         :param even_edges: True for even edges, False for odd edges.
#         :type even_edges: bool
#         """
#         self.n_loops = n_loops
#         self.n_marked_edges = n_marked_edges
#         self.n_hairs = n_hairs
#         self.even_edges = even_edges
#         self.sub_type = sub_types.get(even_edges)
#         max_vertices = 2*n_loops-2 + n_hairs
#         min_vertices = n_marked_edges+1
#         super(ForestedDegSlice, self).__init__(
#             [ForestedGVS(v, n_loops, n_marked_edges, n_hairs, even_edges)
#              for v in range(min_vertices, max_vertices + 1)],
#             n_marked_edges)

#     def get_ordered_param_dict(self):
#         return Shared.OrderedDict([('loops', self.n_loops), ('marked_edges', self.n_marked_edges), ('hairs', self.n_hairs)])

#     def __eq__(self, other):
#         return self.n_loops == other.n_loops \
#             and self.n_marked_edges == other.n_marked_edges and self.n_hairs == other.n_hairs

#     def __str__(self):
#         return ("ForestedDegSlice_%s_%s_%s" % self.get_ordered_param_dict().get_value_tuple()) + self.sub_type

#     def __hash__(self):
#         return hash(str(self))

#     def get_info_plot_path(self):
#         s = "info_vertex_loop_degree_slice_deg_%d_%d_%d_%s_%s" % (
#             self.n_loops, self.n_marked_edges, self.n_hairs, graph_type, self.sub_type)
#         return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)

#     def get_n(self):
#         return self.n_hairs

#     def get_isotypical_projector(self, rep_index):
#         """Returns the SymmetricProjectionOperator corresponding to the isotypical component corresponding to
#         the rep_index-th irrep (as in Partitions(n))
#         """
#         return SymmProjectorDegSlice(self, rep_index)


# class SymmProjectorDegSlice(SymmetricGraphComplex.SymmetricProjectionOperatorDegSlice):
#     def __init__(self, domain, rep_index):
#         """
#         : param domain: Domain vector space of the operator.
#         : type domain: HairyGraphVS
#         : param rep_index: The index of the representation in the list produced by Partitions(h).
#         : type rep_index: int
#         """
#         self.sub_type = domain.sub_type

#         super(SymmProjectorDegSlice, self).__init__(domain, rep_index)

#     def get_ordered_param_dict2(self):
#         do = self.domain
#         return Shared.OrderedDict([('loops', do.n_loops), ('marked edges', do.n_marked_edges), ('hairs', do.n_hairs), ('rep_index', self.rep_index)])

#     def get_matrix_file_path(self):
#         s = "projectionODegSlice%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
#         return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

#     def get_rank_file_path(self):
#         s = "projectionODegSlice%d_%d_%d_%d_rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
#         return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

#     def get_ref_matrix_file_path(self):
#         s = "projectionODegSlice%d_%d_%d_%d.txt" % self.get_ordered_param_dict2().get_value_tuple()
#         return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

#     def get_ref_rank_file_path(self):
#         s = "projectionODegSlice%d_%d_%d_%d.txt.rank.txt" % self.get_ordered_param_dict2().get_value_tuple()
#         return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)


# class ContractEpsToOmegaBiOM(SymmetricGraphComplex.SymmetricBiOperatorMatrix):
#     """Bi operator matrix based on the differentials contract edges and unmark edges.

#     Attributes:
#             - sub_type (str): Sub type of graphs.
#     """

#     def __init__(self, domain, target):
#         self.sub_type = domain.sub_type
#         super(ContractEpsToOmegaBiOM, self).__init__(domain, target, ContractEdgesGO,
#                                                  EpsToOmegaGO)

#     @classmethod
#     def generate_operator(cls, n_loops, n_marked_edges, n_hairs, even_edges):
#         domain = ForestedDegSlice(n_loops,
#                                   n_marked_edges, n_hairs, even_edges)
#         target = ForestedDegSlice(n_loops,
#                                   n_marked_edges-1, n_hairs, even_edges)
#         return cls(domain, target)

#     @staticmethod
#     def is_match(domain, target):
#         """Check whether domain and target degree slices match to generate a corresponding bi operator matrix.

#         The bi operator reduces the degree by one.

#         :param domain: Potential domain vector space of the operator.
#         :type domain: ForestedDegSlice
#         :param target: Potential target vector space of the operator.
#         :type target: ForestedDegSlice
#         :return: bool: True if domain and target match to generate a corresponding bi operator matrix.
#         :rtype: bool
#         """
#         return domain.n_marked_edges - 1 == target.n_marked_edges \
#             and domain.n_loops == target.n_loops \
#             and domain.n_hairs == target.n_hairs \
#             and domain.even_edges == target.even_edges

#     def get_matrix_file_path(self):
#         s = "bi_D_contract_unmark_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
#         return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

#     def get_rank_file_path(self):
#         s = "bi_D_contract_unmark_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
#         return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

#     def restrict_to_isotypical_component(self, rep_index):
#         return RestrictedContractUnmarkGO(self, rep_index)

# class ContractEpsToOmegaD(GraphOperator.Differential):
#     """Differential on the bi graded vector space based on the operators contract edges and delete edges.

#     Only for graphs with odd edges.
#     """

#     def __init__(self, graded_sum_vs):
#         """Initialize the contract and delete edges differential with the underlying bi graded vector space.

#         :param graded_sum_vs: Underlying bi graded vector space.
#         :type graded_sum_vs: VertexLoopBigradedSumVS
#         """
#         super(ContractEpsToOmegaD, self).__init__(graded_sum_vs,
#                                               ContractEpsToOmegaBiOM.generate_op_matrix_list(graded_sum_vs))

#     def get_type(self):
#         return 'contract edges and unmark edges'

#     def get_cohomology_plot_path(self):
#         sub_type = self.sum_vector_space.sub_type
#         s = "cohomology_dim_contract_edges_unmark_edges_D_%s_%s" % (
#             graph_type, sub_type)
#         return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

#     def get_info_plot_path(self):
#         sub_type = self.sum_vector_space.sub_type
#         s = "info_contract_edges_unmark_edges_D_%s_%s" % (
#             graph_type, sub_type)
#         return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

#     def get_ordered_cohomology_param_range_dict(self):
#         s = self.sum_vector_space
#         return Shared.OrderedDict([('loops', s.l_range), ('marked_edges', s.m_range), ('hairs', s.h_range)])


# class RestrictedContractUnmarkD(SymmetricGraphComplex.SymmetricDifferential):
#     def get_type(self):
#         return 'isotypical contract edges'

#     def get_cohomology_plot_path(self):
#         sub_type = self.diff.sum_vector_space.sub_type
#         s = "cohomology_dim_contract_edges_unmark_edges_D_iso_%s_%s" % (
#             graph_type, sub_type)
#         return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

#     def get_info_plot_path(self):
#         sub_type = self.diff.sum_vector_space.sub_type
#         s = "info_contract_edges_unmark_edges_D_iso_%s_%s" % (
#             graph_type, sub_type)
#         return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


# class ForestedContractUnmarkBiGC(GraphComplex.GraphComplex):
#     """Bi complex based on ordinary simple graphs and the differentials contract edges and unmark edges.

#     Attributes:
#         - deg_range (range): Range for the total degree.
#         - even_edges (bool): True for even edges, False for odd edges.
#         - sub_type (str): Sub type of graphs.
#     """

#     def __init__(self, l_range, m_range, h_range, even_edges, isotypical=False):
#         """Initialize the bi complex.

#         :param deg_range: Range for the degree.
#         :type deg_range: range
#         :param even_edges: True for even edges, False for odd edges.
#         :type even_edges: bool
#         """
#         self.l_range = l_range
#         self.m_range = m_range
#         self.h_range = h_range
#         self.even_edges = even_edges
#         self.sub_type = sub_types.get(self.even_edges)
#         graded_sum_vs = ForestedBigradedSumVS(
#             l_range, m_range, h_range, self.even_edges)
#         self.contract_unmarkD = ContractUnmarkD(graded_sum_vs)
#         if isotypical:
#             super(ForestedContractUnmarkBiGC, self).__init__(
#                 graded_sum_vs, [RestrictedContractUnmarkD(self.contract_unmarkD)])
#         else:
#             super(ForestedContractUnmarkBiGC, self).__init__(
#                 graded_sum_vs, [self.contract_unmarkD])

#     def __str__(self):
#         return '<%s graphs bi-complex with %s>' % (graph_type, str(self.sub_type))

#     def print_dim_and_eulerchar(self):
#         for h in self.h_range:
#             for l in self.l_range:
#                 ds = [ForestedDegSlice(l, m, h, self.even_edges).get_dimension()
#                       for m in self.m_range]
#                 eul = sum([(1 if j % 2 == 0 else -1) *
#                            d for j, d in enumerate(ds)])
#                 print("Dimensions (h,l) ",
#                       h, l, self.sub_type, ":", ds, "Euler", eul)

#     def print_cohomology_dim(self):
#         for h in self.h_range:
#             for l in self.l_range:
#                 cohomdict = {}
#                 for m in self.m_range:
#                     D1 = ContractUnmarkBiOM.generate_operator(
#                         l, m, h, self.even_edges)
#                     D2 = ContractUnmarkBiOM.generate_operator(
#                         l, m+1, h, self.even_edges)
#                     try:
#                         d = ForestedDegSlice(
#                             l, m, h, self.even_edges).get_dimension()
#                         r1 = D1.get_matrix_rank()
#                         r2 = D2.get_matrix_rank()
#                         cohomdict[m] = d-r1-r2
#                     except:
#                         pass

#                 print("Cohomology Dimensions (h,l) ",
#                       h, l, self.sub_type, ":", cohomdict)

#     def build_basis(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
#         # print("Building auxiliary pregraphs...")
#         # PreForestedGraphSumVS.compute_all_pregraphs(-1,
#         #                                             max(self.l_range), max(self.m_range), max(self.h_range), self.even_edges, ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar, info_tracker=info_tracker)
#         # print("Done.")
#         return super().build_basis(ignore_existing_files, n_jobs, progress_bar, info_tracker)
