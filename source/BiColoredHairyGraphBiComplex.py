"""Simple bicolored hairy graph bi complex based on the differentials contract edges and split edges.
Only for graphs with odd edges and even hairs for both hair colors."""

__all__ = ['ContractSplitBiOM', 'VertexLoopDegSlice', 'VertexLoopBigradedSumVS', 'ContractSplitD',
           'BiColoredHairyContractSplitBiGC']

import os
import itertools
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared as Shared
import Parameters
import BiColoredHairyGraphComplex


class ContractSplitBiOM(GraphOperator.BiOperatorMatrix):
    def __init__(self, domain, target):
        """Bi operator matrix based on the differentials contract edges and split edges.

        Attributes:
                - sub_type (str): Sub type of graphs.
        """
        self.sub_type = domain.sub_type
        super(ContractSplitBiOM, self).__init__(domain, target, BiColoredHairyGraphComplex.ContractEdgesGO,
                                                BiColoredHairyGraphComplex.SplitEdgesGO)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target degree slices match to generate a corresponding bi operator matrix.

        The bi operator reduces the degree by one and increases the minimal number of hairs by one for both hair colors.

        :param domain: Potential domain vector space of the operator.
        :type domain: VertexLoopDegSlice
        :param target: Potential target vector space of the operator.
        :type target: VertexLoopDegSlice
        :return: True if domain and target match to generate a corresponding bi operator matrix.
        :rtype: bool
        """
        return domain.deg - 1 == target.deg and domain.h_a_min + 1 == target.h_a_min and \
            domain.h_b_min + 1 == target.h_b_min

    def get_matrix_file_path(self):
        s = "bi_D_contract_split_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, BiColoredHairyGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_contract_split_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, BiColoredHairyGraphComplex.graph_type, self.sub_type, s)


class VertexLoopDegSlice(GraphVectorSpace.DegSlice):
    """Degree slice of hairy graphs with the two degrees number of vertices and loops.

    Total degree = n_vertices + n_loops

    Attributes:
        - h_a_min (int): Minimal number of hairs_a. Can be negative.
        - h_b_min (int): Minimal number of hairs_b. Can be negative.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs_a (bool): True for even hairs_a, False for odd hairs_a.
        - even_hairs_b (bool): True for even hairs_a, False for odd hairs_b.
    """

    def __init__(self, deg, h_a_min, h_b_min, even_edges, even_hairs_a, even_hairs_b):
        """Initialize the degree slice.

        :param deg: Total degree of the degree slice.
        :type deg: int
        :param h_a_min: Minimal number of hairs_a.
        :type h_a_min: int
        :param h_b_min: Minimal number of hairs_b.
        :type h_b_min: int
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: True for even hairs_a, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: True for even hairs_b, False for odd hairs_b.
        :type even_hairs_b: bool
        """
        self.h_a_min = h_a_min
        self.h_b_min = h_b_min
        self.sub_type = BiColoredHairyGraphComplex.get_sub_type(
            even_edges, h_a_min, h_b_min)
        super(VertexLoopDegSlice, self).__init__(
            [BiColoredHairyGraphComplex.BiColoredHairyGraphVS(v, deg - v, self.h_a_min + v, self.h_b_min + v,
                                                              even_edges, even_hairs_a, even_hairs_b)
             for v in range(0, deg + 1)], deg)

    def __hash__(self):
        return hash(str(self))

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('deg', self.deg), ('min_hairs_a', self.h_a_min), ('min_hairs_b', self.h_b_min)])

    def __eq__(self, other):
        return self.deg == other.deg and self.h_a_min == other.h_a_min and self.h_b_min == other.h_b_min

    def get_info_plot_path(self):
        s = "info_vertex_loop_degree_slice_deg_%d_h_a_min_%d_h_b_min_%d_%s_%s" % (self.deg, self.h_a_min, self.h_b_min,
                                                                                  BiColoredHairyGraphComplex.graph_type,
                                                                                  self.sub_type)
        return os.path.join(Parameters.plots_dir, BiColoredHairyGraphComplex.graph_type, self.sub_type, s)


class VertexLoopBigradedSumVS(GraphVectorSpace.SumVectorSpace):
    """Bi graded vector space based on simple hairy graphs with two colors of hairs.

    Bi grading according to the number of vertices and loops.
    Direct sum of degree slices.

    Attributes:
        - deg_range (range): Range for the total degree.
        - h_a_min_range (range): Range for minimal number of hairs_a for the degree slice with the highest degree.
        - h_b_min_range (range): Range for minimal number of hairs_b for the degree slice with the highest degree.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs_a (bool): True for even hairs_a, False for odd hairs_a.
        - even_hairs_b (bool): True for even hairs_b, False for odd hairs_b.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, deg_range, h_a_min_range, h_b_min_range, even_edges, even_hairs_a, even_hairs_b):
        """ Initialize the bi graded vector space.

        :param deg_range: Range for the degree.
        :type deg_range: range
        :param h_a_min_range: Range for minimal number of hairs_a for the degree slice with the highest degree.
        :type h_a_min_range: range
        :param h_b_min_range: Range for minimal number of hairs_b for the degree slice with the highest degree.
        :type h_b_min_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: True for even hairs_a, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: True for even hairs_b, False for odd hairs_b.
        :type even_hairs_b: bool
        """
        self.deg_range = deg_range
        self.h_a_min_range = h_a_min_range
        self.h_b_min_range = h_b_min_range
        self.sub_type = BiColoredHairyGraphComplex.get_sub_type(
            even_edges, even_hairs_a, even_hairs_b)
        max_deg = max(self.deg_range)

        vs_list = []
        for (deg, h_a_min, h_b_min) in itertools.product(self.deg_range, self.h_a_min_range, self.h_b_min_range):
            if even_hairs_a == even_hairs_b and h_a_min < h_b_min:
                continue  # Symmetry between a and b hairs.
            vs_list.append(VertexLoopDegSlice(deg, h_a_min + (max_deg - deg), h_b_min + (max_deg - deg), even_edges,
                                              even_hairs_a, even_hairs_b))
        super(VertexLoopBigradedSumVS, self).__init__(vs_list)

    def get_type(self):
        return '%s graphs with %s' % (BiColoredHairyGraphComplex.graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('deg', self.deg_range), ('min_hairs_a', self.h_a_min_range),
                                   ('min_hairs_b', self.h_b_min_range)])

    def get_info_plot_path(self):
        s = "info_vertex_loop_bigraded_vector_space_%s_%s" % (
            BiColoredHairyGraphComplex.graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, BiColoredHairyGraphComplex.graph_type, self.sub_type, s)


class ContractSplitD(GraphOperator.Differential):
    """Differential on the bi graded vector space based on the operators contract edges and split edges.

    Only for graphs with odd edges and even hairs for both hair colors.
    """

    def __init__(self, graded_sum_vs):
        """Initialize the contract and split edges differential with the underlying bi graded vector space.

        :param graded_sum_vs: Underlying bi graded vector space.
        :type graded_sum_vs: VertexLoopBigradedSumVS
        """
        super(ContractSplitD, self).__init__(graded_sum_vs,
                                             ContractSplitBiOM.generate_op_matrix_list(graded_sum_vs))

    def get_type(self):
        return 'contract edges and split edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_split_edges_D_%s_%s" % (
            BiColoredHairyGraphComplex.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, BiColoredHairyGraphComplex.graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_split_edges_D_%s_%s" % (
            BiColoredHairyGraphComplex.graph_type, sub_type)
        return os.path.join(Parameters.web_dir, BiColoredHairyGraphComplex.graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_edges_split_edges_D_%s_%s" % (
            BiColoredHairyGraphComplex.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, BiColoredHairyGraphComplex.graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        deg_range = self.sum_vector_space.deg_range
        h_a_min_min = min(self.sum_vector_space.h_a_min_range)
        h_a_min_max = max(self.sum_vector_space.h_a_min_range) + \
            (max(deg_range) - min(deg_range))
        h_a_min_range = range(h_a_min_min, h_a_min_max + 1)
        h_b_min_min = min(self.sum_vector_space.h_a_min_range)
        h_b_min_max = max(self.sum_vector_space.h_a_min_range) + \
            (max(deg_range) - min(deg_range))
        h_b_min_range = range(h_b_min_min, h_b_min_max + 1)
        return Shared.OrderedDict([('deg', deg_range), ('min_hairs_a', h_a_min_range),
                                   ('min_hairs_b', h_b_min_range)])

    def get_cohomology_plot_parameter_order(self):
        return (1, 2, 0)


class BiColoredHairyContractSplitBiGC(GraphComplex.GraphComplex):
    """Bi complex based on simple hairy graphs and the differentials contract edges and split edges.

    Only for graphs with odd edges and even hairs for both hair colors.

    Attributes:
        - deg_range (range): Range for the total degree.
        - h_a_min_range (range): Range for minimal number of hairs_a for the degree slice with the highest degree.
        - h_b_min_range (range): Range for minimal number of hairs_b for the degree slice with the highest degree.
        - even_edges (bool): True for even edges, False for odd edges.
        - even_hairs_a (bool): True for even hairs_a, False for odd hairs_a.
        - even_hairs_b (bool): True for even hairs_b, False for odd hairs_b.
        - sub_type (str): Sub type of graphs.
    """

    def __init__(self, deg_range, h_a_min_range, h_b_min_range, even_edges, even_hairs_a, even_hairs_b):
        """Initialize the bi complex.

        :param deg_range: Range for the degree.
        :type deg_range: range
        :param h_a_min_range: Range for minimal number of hairs_a for the degree slice with the highest degree.
        :type h_a_min_range: range
        :param h_b_min_range: Range for minimal number of hairs_b for the degree slice with the highest degree.
        :type h_b_min_range: range
        :param even_edges: True for even edges, False for odd edges.
        :type even_edges: bool
        :param even_hairs_a: True for even hairs_a, False for odd hairs_a.
        :type even_hairs_a: bool
        :param even_hairs_b: True for even hairs_b, False for odd hairs_b.
        :type even_hairs_b: bool
        """
        self.sub_type = BiColoredHairyGraphComplex.get_sub_type(
            even_edges, even_hairs_a, even_hairs_b)

        graded_sum_vs = VertexLoopBigradedSumVS(deg_range, h_a_min_range, h_b_min_range, even_edges, even_hairs_a,
                                                even_hairs_b)
        super(BiColoredHairyContractSplitBiGC, self).__init__(
            graded_sum_vs, [ContractSplitD(graded_sum_vs)])

    def __str__(self):
        return '<%s graphs bi-complex with %s>' % (BiColoredHairyGraphComplex.graph_type, str(self.sub_type))
