"""Graph complexes based on ordinary (simple, without multi edges) graphs.
Implemented Differentials: Contract edges, delete edges."""


__all__ = ['graph_type', 'sub_types', 'OrdinaryGVS', 'OrdinaryGraphSumVS', 'ContractEdgesGO', 'ContractEdgesD',
           'DeleteEdgesGO', 'DeleteEdgesD', 'OrdinaryGC']

import itertools
from sage.all import *
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import NautyInterface
import Parameters


graph_type = "ordinary"

sub_types = {True: "even_edges", False: "odd_edges"}


# ------- Graph Vector Space --------
class OrdinaryGVS(GraphVectorSpace.GraphVectorSpace):
    """Ordinary graph vector space.

    Sub vector space with specified number of vertices, loops and even or odd edges and at least trivalent vertices.
    No multiple edges.

    Attributes:
        n_vertices (non-negative int): Number of vertices.

        n_loops (non-negative int): Number of loops.

        even_edges (bool): True for even edges, False for odd edges.

        n_edges (non-negative int): Number of edges.

        sub_type (str): Sub type of graphs.

    """
    def __init__(self, n_vertices, n_loops, even_edges):
        """Initialize the ordinary graph vector space.

        :param n_vertices: non-negative int: Number of vertices.
        :param n_loops: non-negative int: Number of loops.
        :param even_edges: bool: True for even edges, False for odd edges.
        """
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        self.sub_type = sub_types.get(self.even_edges)
        super(OrdinaryGVS, self).__init__()

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def __eq__(self, other):
        return self.n_vertices == other.n_vertices and self.n_loops == other.n_loops \
               and self.even_edges == other.even_edges

    def get_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_basis_file_path(self):
        s = "gra%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])

    def get_partition(self):
        return None

    def is_valid(self):
        # Vertices at least trivalent. Positive number of vertices. Non-negative number of loops.
        # At most fully connected graph, no multiple edges.
        return (3 * self.n_vertices <= 2 * self.n_edges) and self.n_vertices > 0 and self.n_loops >= 0 \
               and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2

    def get_work_estimate(self):
        # Returns the number of possible graphs as work estimate.
        if not self.is_valid():
            return 0
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def get_generating_graphs(self):
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            return []
        return NautyInterface.list_simple_graphs(self.n_vertices, self.n_edges)

    def perm_sign(self, G, p):
        if self.even_edges:
            # The sign is (induced sign on vertices) * (induced sign edge orientations)
            sign = Shared.Perm(p).signature()
            for (u, v) in G.edges(labels=False):
                # We assume the edge is always directed from the larger to smaller index
                if (u < v and p[u] > p[v]) or (u > v and p[u] < p[v]):
                    sign *= -1
            return sign
        else:
            # The sign is (induced sign of the edge permutation)
            # We assume the edges are always lexicographically ordered
            # For the computation we use that G.edges() returns the edges in lex ordering
            # We first label the edges on a copy of G lexicographically
            G1 = copy(G)
            Shared.enumerate_edges(G1)
            # We permute the graph, and read of the new labels
            G1.relabel(p, inplace=True)
            return Shared.Perm([j for (u, v, j) in G1.edges()]).signature()


class OrdinaryGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of ordinary graph vector spaces with specified edge parity.

    Attributes:
        v_range (range): Range for the number of vertices.

        l_range (range): Range for the number of loops.

        even_edges (bool): True for even edges, False for odd edges.

        sub_type (str): Sub type of graphs.
    """
    def __init__(self, v_range, l_range, even_edges):
        """Initialize the sum vector space.

        :param v_range: range: Range for the number of vertices.
        :param l_range: range: Range for the number of loops.
        :param even_edges: bool: True for even edges, False for odd edges.
        """
        self.v_range = v_range
        self.l_range = l_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        vs_list = [OrdinaryGVS(v, l, self.even_edges) for (v, l) in itertools.product(self.v_range, self.l_range)]
        super(OrdinaryGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range)])


# ------- Operators --------
class ContractEdgesGO(GraphOperator.GraphOperator):
    """Contract edges graph operator.

    Operates on an ordinary graph by contracting an edge and unifying the two adjacent vertices.

    Attributes:
        sub_type (str): Graphs sub type of the domain.
    """
    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the contract edges graph operator.

        :param domain: OrdinaryGVS: Domain vector space of the operator.
        :param target: OrdinaryGVS: Target vector space of the operator.
        """
        if not ContractEdgesGO.is_match(domain, target):
            raise ValueError("Domain and target not consistent for contract edges operator")
        self.sub_type = domain.sub_type
        super(ContractEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding contract edges graph operator.

        :param domain: OrdinaryGVS: Potential domain vector space of the operator.
        :param target: OrdinaryGVS: Potential target vector space of the operator.
        :return: bool: True if domain and target match to generate a corresponding contract edges graph operator.
        """
        return domain.n_vertices - 1 == target.n_vertices and domain.n_loops == target.n_loops \
                and domain.even_edges == target.even_edges

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, even_edges):
        """Returns a contract edge graph operator.

        :param n_vertices: non-negative int: Number of vertices of the domain.
        :param n_loops: non-negative int: Number of loops of the domain.
        :param even_edges: bool: True for even edges, False for odd edges.
        :return: ContractEdgesGO: Contract edges graph operator based on the specified domain vector space.
        """
        domain = OrdinaryGVS(n_vertices, n_loops, even_edges)
        target = OrdinaryGVS(n_vertices - 1, n_loops, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "contractD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        # Operates on the graph G by contracting an edge and unifying the adjacent vertices.
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            pp = Shared.permute_to_left((u, v), range(0, self.domain.n_vertices))
            sgn = self.domain.perm_sign(G, pp)
            G1 = copy(G)
            G1.relabel(pp, inplace=True)
            Shared.enumerate_edges(G1)
            previous_size = G1.size()
            G1.merge_vertices([0,1])
            if (previous_size - G1.size()) != 1:
                continue
            G1.relabel(list(range(0,G1.order())), inplace=True)
            if not self.domain.even_edges:
                p = [j for (a, b, j) in G1.edges()]
                sgn *= Permutation(p).signature()
            else:
                sgn *= -1 #TODO overall sign for even edges
            image.append((G1, sgn))
        return image


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""
    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space.

        :param sum_vector_space: OrdinaryGraphSumVS: Underlying vector space.
        """
        super(ContractEdgesD, self).__init__(sum_vector_space, ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.get_vs_list()[0].sub_type
        s = "cohomology_dim_contrct_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


class DeleteEdgesGO(GraphOperator.GraphOperator):
    """Delete edges graph operator.

    Operates on an ordinary graph by deleting an edge.
    Only for graphs with odd edges.

    Attributes:
        sub_type (str): Graphs sub type of the domain.
    """
    def __init__(self, domain, target):
        """Initialize the domain and target vector space of the delete edges graph operator.

        :param domain: OrdinaryGVS: Domain vector space of the operator.
        :param target: OrdinaryGVS: Target vector space of the operator.
        """
        if not DeleteEdgesGO.is_match(domain, target):
            raise ValueError("Domain and target not consistent for delete edges operator")
        self.sub_type = domain.sub_type
        super(DeleteEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding delete edges graph operator.

        :param domain: OrdinaryGVS: Potential domain vector space of the operator.
        :param target: OrdinaryGVS: Potential target vector space of the operator.
        :return: bool: True if domain and target match to generate a corresponding contract edges graph operator.
        """
        return domain.n_vertices == target.n_vertices and domain.n_loops - 1 == target.n_loops \
                and domain.even_edges == target.even_edges

    @classmethod
    def generate_operator(cls, n_vertices, n_loops, even_edges):
        """Returns a delete edge graph operator.

        :param n_vertices: non-negative int: Number of vertices of the domain.
        :param n_loops: non-negative int: Number of loops of the domain.
        :param even_edges: bool: True for even edges, False for odd edges.
        :return: ContractEdgesGO: Contract edges graph operator based on the specified domain vector space.
        """
        domain = OrdinaryGVS(n_vertices, n_loops, even_edges)
        target = OrdinaryGVS(n_vertices, n_loops - 1, even_edges)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "delete_edgesD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "delete_edgesD%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        s = "delete_edgesD%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, self.sub_type, s)

    def get_ref_rank_file_path(self):
        s = "delete_edgesD%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
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
        return 'delete edges'

    def operate_on(self, G):
        # Operates on the graph G by deleting an edge.
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            G1 = copy(G)
            G1.delete_edge((u, v))
            sgn = -1 if i % 2 else 1
            image.append((G1, sgn))
        return image


class DeleteEdgesD(GraphOperator.Differential):
    """Delete edges differential.

    Only for graphs with odd edges.
    """
    def __init__(self, sum_vector_space):
        """Initialize the delete edges differential with the underlying sum vector space.

        :param sum_vector_space: OrdinaryGraphSumVS: Underlying vector space.
        """
        super(DeleteEdgesD, self).__init__(sum_vector_space, DeleteEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'delete edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.get_vs_list()[0].sub_type
        s = "cohomology_dim_delete_edges_D_%s_%s" % (graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)


# ------- Graph Complexes --------
class OrdinaryGC(GraphComplex.GraphComplex):
    """Graph complex for ordinary graphs.

    Attributes:
        v_range (range): Range for the number of vertices.

        l_range (range): Range for the number of loops.

        even_edges (bool): True for even edges, False for odd edges.

        sub_type (str): Sub type of graphs.
    """
    def __init__(self, v_range, l_range, even_edges, differentials):
        """Initialize the graph complex.

        :param v_range: range: Range for the number of vertices.
        :param l_range: range: Range for the number of loops.
        :param even_edges: bool: True for even edges, False for odd edges.
        :param differentials: list(str): List of differentials. Options: 'contract', 'delete'.
        """
        self.v_range = v_range
        self.l_range = l_range
        self.even_edges = even_edges
        self.sub_type = sub_types.get(self.even_edges)

        sum_vector_space = OrdinaryGraphSumVS(v_range, l_range, even_edges)
        differential_list = []
        if not differentials <= {'contract', 'delete'}:
            raise ValueError("Differentials for ordinary graph complex: 'contract', 'delete'")
        if 'contract' in differentials:
            contract_edges_dif = ContractEdgesD(sum_vector_space)
            differential_list.append(contract_edges_dif)
        if 'delete' in differentials:
            delete_edges_dif = DeleteEdgesD(sum_vector_space)
            differential_list.append(delete_edges_dif)
        super(OrdinaryGC, self).__init__(sum_vector_space, differential_list)

    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))
