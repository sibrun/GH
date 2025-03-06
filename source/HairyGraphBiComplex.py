"""Simple hairy graph bi complex based on the differentials contract edges and edge to one hair.
Only for graphs with odd edges and odd hairs."""

___all__ = ['CeEt1hBiOM', 'VertexLoopDegSlice', 'VertexLoopBigradedSumVS', 'CeEt1hD', 'HairyCeEt1hBiGC' ]

import os
import itertools
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared as Shared
import Parameters
import HairyGraphComplex


class CeEt1hBiOM(GraphOperator.BiOperatorMatrix):
    """Bi operator matrix based on the differentials contract edges and edge to one hair.

    Attributes:
            - sub_type (str): Sub type of graphs.
    """
    def __init__(self, domain, target):
        self.sub_type = domain.sub_type
        super().__init__(domain, target, HairyGraphComplex.ContractEdgesGO,
                                         HairyGraphComplex.EdgeToOneHairGO)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target degree slices match to generate a corresponding bi operator matrix.

        The bi operator reduces the degree by one and increases the minimal number of hairs by one.

        :param domain: Potential domain vector space of the operator.
        :type domain: VertexLoopDegSlice
        :param target: Potential target vector space of the operator.
        :type target: VertexLoopDegSlice
        :return: True if domain and target match to generate a corresponding bi operator matrix.
        :rtype: bool
        """
        return domain.deg - 1 == target.deg and domain.h_min + 1 == target.h_min

    def get_matrix_file_path(self):
        s = "bi_D_ce_et1h_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, HairyGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_ce_et1h_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, HairyGraphComplex.graph_type, self.sub_type, s)


class VertexLoopDegSlice(GraphVectorSpace.DegSlice):
    """Degree slice of hairy graphs with the two degrees number of vertices and loops.

    Total degree = n_vertices + n_loops

    Attributes:
        - h_min (int): Minimal number of hairs. Can be negative.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs (bool): True for even hairs, False for odd hairs.
    """
    def __init__(self, deg, h_min, even_edges, even_hairs):
        """Initialize the degree slice.

        :param deg: Total degree of the degree slice.
        :type deg: int
        :param h_min: Minimal number of hairs.
        :type h_min: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs: True for even hairs, False for odd hairs.
        :type even_hairs: bool
        """
        self.h_min = h_min
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = HairyGraphComplex.sub_types.get((self.even_edges, self.even_hairs))
        super().__init__(
            [HairyGraphComplex.HairyGraphVS(v, deg - v, self.h_min + v, self.even_edges, self.even_hairs)
             for v in range(0, deg + 1)], deg)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('deg', self.deg), ('min_hairs', self.h_min)])

    def __eq__(self, other):
        return self.deg == other.deg and self.h_min == other.h_min

    def get_info_plot_path(self):
        s = "info_vertex_loop_degree_slice_deg_%d_h_min_%d_%s_%s" % (self.deg, self.h_min, HairyGraphComplex.graph_type,
                                                                     self.sub_type)
        return os.path.join(Parameters.plots_dir, HairyGraphComplex.graph_type, self.sub_type, s)


class VertexLoopBigradedSumVS(GraphVectorSpace.SumVectorSpace):
    """Bi graded vector space based on simple hairy graphs.

    Bi grading according to the number of vertices and loops.
    Direct sum of degree slices.

    Attributes:
        - deg_range (range): Range for the total degree.
        - h_min_range (range): Range for minimal number of hairs for the degree slice with the highest degree.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs (bool): True for even hairs, False for odd hairs.
        - sub_type (str): Sub type of graphs.
    """
    def __init__(self, deg_range, h_min_range, even_edges, even_hairs):
        """ Initialize the bi graded vector space.

        :param deg_range: Range for the degree.
        :type deg_range: range
        :param h_min_range: Range for minimal number of hairs for the degree slice with the highest degree.
        :type h_min_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs: True for even hairs, False for odd hairs.
        :type even_hairs: bool
        """
        self.deg_range = deg_range
        self.h_min_range = h_min_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = HairyGraphComplex.sub_types.get((self.even_edges, self.even_hairs))
        max_deg = max(self.deg_range)
        super().__init__(
            [VertexLoopDegSlice(deg, h_min + (max_deg - deg), self.even_edges, self.even_hairs) for (deg, h_min) in
             itertools.product(self.deg_range, self.h_min_range)])

    def get_type(self):
        return '%s graphs with %s' % (HairyGraphComplex.graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('deg', self.deg_range), ('min_hairs', self.h_min_range)])

    def get_info_plot_path(self):
        s = "info_vertex_loop_bigraded_vector_space_%s_%s" % (HairyGraphComplex.graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, HairyGraphComplex.graph_type, self.sub_type, s)


class ContractEdgeToOneHD(GraphOperator.Differential):
    """Differential on the bi graded vector space based on the operators contract edges and edge to one hair.

    Only for graphs with odd edges and odd hairs.
    """
    def __init__(self, graded_sum_vs):
        """Initialize the contract and edge to one hair differential with the underlying bi graded vector space.

        :param graded_sum_vs: Underlying bi graded vector space.
        :type graded_sum_vs: VertexLoopBigradedSumVS
        """
        super().__init__(graded_sum_vs, CeEt1hBiOM.generate_op_matrix_list(graded_sum_vs))

    def get_type(self):
        return 'contract edges and edge to one hair'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_edge_to_one_hair_D_%s_%s" % (HairyGraphComplex.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, HairyGraphComplex.graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_edge_to_one_hair_D_%s_%s" % (HairyGraphComplex.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, HairyGraphComplex.graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        deg_range = self.sum_vector_space.deg_range
        h_min_min = min(self.sum_vector_space.h_min_range)
        h_min_max = max(self.sum_vector_space.h_min_range) + (max(deg_range) - min(deg_range))
        h_range = range(h_min_min, h_min_max + 1)
        return Shared.OrderedDict([('deg', deg_range), ('min_hairs', h_range)])

    def get_cohomology_plot_parameter_order(self):
        return (0, 1)


class HairyCeEt1hBiGC(GraphComplex.GraphComplex):
    """Bi complex based on simple hairy graphs and the differentials contract edges and edge to one hair.

    Only for graphs with odd edges and odd hairs.

    Attributes:
        - deg_range (range): Range for the total degree.
        - h_min_range (range): Range for minimal number of hairs for the degree slice with the highest degree.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs (bool): True for even hairs, False for odd hairs.
        - sub_type (str): Sub type of graphs.
    """
    def __init__(self, deg_range, h_min_range, even_edges, even_hairs):
        """Initialize the bi complex.

        :param deg_range: Range for the degree.
        :type deg_range: range
        :param h_min_range: Range for minimal number of hairs for the degree slice with the highest degree.
        :type h_min_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs: True for even hairs, False for odd hairs.
        :type even_hairs: bool
        """
        self.deg_range = deg_range
        self.h_min_range = h_min_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = HairyGraphComplex.sub_types.get((self.even_edges, self.even_hairs))

        graded_sum_vs = VertexLoopBigradedSumVS(self.deg_range, self.h_min_range, self.even_edges, self.even_hairs)
        super().__init__(graded_sum_vs, [ContractEdgeToOneHD(graded_sum_vs)])

    def __str__(self):
        return '<%s graphs bi-complex with %s>' % (HairyGraphComplex.graph_type, str(self.sub_type))

