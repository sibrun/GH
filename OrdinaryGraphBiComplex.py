"""Ordinary simple graph bi complex based on the differentials contract edges and delete edges.
Only for graphs with odd edges."""

__all__ = ['ContractDeleteBiOM', 'VertexLoopDegSlice', 'VertexLoopBigradedSumVS', 'ContractDeleteD',
           'OrdinaryContractDeleteBiGC']

import os
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import Parameters
import OrdinaryGraphComplex


class ContractDeleteBiOM(GraphOperator.BiOperatorMatrix):
    """Bi operator matrix based on the differentials contract edges and delete edges.

    Attributes:
            sub_type (str): Sub type of graphs.
    """
    def __init__(self, domain, target):
        self.sub_type = domain.get_vs_list()[0].sub_type
        super(ContractDeleteBiOM, self).__init__(domain, target, OrdinaryGraphComplex.ContractEdgesGO,
                                                 OrdinaryGraphComplex.DeleteEdgesGO)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target degree slices match to generate a corresponding bi operator matrix.

        :param domain: VertexLoopDegSlice: Potential domain vector space of the operator.
        :param target: VertexLoopDegSlice: Potential target vector space of the operator.
        :return: bool: True if domain and target match to generate a corresponding bi operator matrix.
        """
        (d_deg,) = domain.get_ordered_param_dict().get_value_tuple()
        (t_deg,) = target.get_ordered_param_dict().get_value_tuple()
        return d_deg - 1 == t_deg

    def get_matrix_file_path(self):
        s = "bi_D_ce_dele_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, OrdinaryGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_ce_dele_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, OrdinaryGraphComplex.graph_type, self.sub_type, s)


class VertexLoopDegSlice(GraphVectorSpace.DegSlice):
    """Degree slice of ordinary graphs with the two degrees number of vertices and loops.

    Degree = n_vertices + n_loops

    Attributes:
        even_edges (bool): True for even edges, False for odd edges.

    """
    def __init__(self, deg, even_edges):
        """Initialize the degree slice.

        :param deg: non-negative int: Degree of the degree slice.
        :param even_edges: bool: True for even edges, False for odd edges.
        """
        self.even_edges = even_edges
        super(VertexLoopDegSlice, self).__init__(
            [OrdinaryGraphComplex.OrdinaryGVS(v, deg - v, self.even_edges) for v in range(0, deg + 1)], deg)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('deg', self.deg)])

    def __eq__(self, other):
        return self.deg == other.deg


class VertexLoopBigradedSumVS(GraphVectorSpace.SumVectorSpace):
    """Bi graded vector space based on ordinary simple graphs.

    Bi grading according to the number of vertices and loops.
    Direct sum of degree slices.

    Attributes:
        deg_range (range): Range for the degree.

        even_edges (bool): True for even edges, False for odd edges.

        sub_type (str): Sub type of graphs.
    """
    def __init__(self, deg_range, even_edges):
        """ Initialize the bi graded vector space.

        :param deg_range: range: Range for the degree.
        :param even_edges: bool: True for even edges, False for odd edges.
        """
        self.deg_range = deg_range
        self.even_edges = even_edges
        self.sub_type = OrdinaryGraphComplex.sub_types.get(even_edges)
        super(VertexLoopBigradedSumVS, self).__init__([VertexLoopDegSlice(deg, self.even_edges) for deg in self.deg_range])

    def get_type(self):
        return '%s graphs with %s' % (OrdinaryGraphComplex.graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('deg', self.deg_range)])


class ContractDeleteD(GraphOperator.Differential):
    """Differential on the bi graded vector space based on the operators contract edges and delete edges.

    Only for odd edges.
    """
    def __init__(self, graded_sum_vs):
        """Initialize the contract and delete edges differential with the underlying bi graded vector space.

        :param graded_sum_vs: VertexLoopBigradedSumVS: Underlying bi graded vector space.
        """
        super(ContractDeleteD, self).__init__(graded_sum_vs, ContractDeleteBiOM.generate_op_matrix_list(graded_sum_vs))

    def get_type(self):
        return 'contract edges and delete edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_delete_edges_D_%s_%s" % (OrdinaryGraphComplex.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, OrdinaryGraphComplex.graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        deg_range = self.sum_vector_space.deg_range
        return Shared.OrderedDict([('deg', deg_range)])


class OrdinaryContractDeleteBiGC(GraphComplex.GraphComplex):
    """Bi complex based on ordinary simple graphs and the differentials contract edges and delete edges.

    Only for odd edges.

    Attributes:
        deg_range (range): Range for the degree.

        even_edges (bool): True for even edges, False for odd edges.

        sub_type (str): Sub type of graphs.
    """
    def __init__(self, deg_range, even_edges):
        """Initialize the bi complex.

        :param deg_range: range: Range for the degree.
        :param even_edges: bool: True for even edges, False for odd edges.
        """
        self.deg_range = deg_range
        self.even_edges = even_edges
        self.sub_type = OrdinaryGraphComplex.sub_types.get(self.even_edges)
        graded_sum_vs = VertexLoopBigradedSumVS(self.deg_range, self.even_edges)
        super(OrdinaryContractDeleteBiGC, self).__init__(graded_sum_vs, [ContractDeleteD(graded_sum_vs)])

    def __str__(self):
        return '<%s graphs bi-complex with %s>' % (OrdinaryGraphComplex.graph_type, str(self.sub_type))