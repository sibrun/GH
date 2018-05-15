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
        self.sub_type = domain.get_vs_list()[0].sub_type
        super(ContractSplitBiOM, self).__init__(domain, target, BiColoredHairyGraphComplex.ContractEdgesGO,
                                                BiColoredHairyGraphComplex.SplitEdgesGO)

    @staticmethod
    def is_match(domain, target):
        (d_deg, d_h_a_min, d_h_b_min) = domain.get_ordered_param_dict().get_value_tuple()
        (t_deg, t_h_a_min, t_h_b_min) = target.get_ordered_param_dict().get_value_tuple()
        return (d_deg - 1, d_h_a_min + 1, d_h_b_min + 1) == (t_deg, t_h_a_min, t_h_b_min)

    def get_matrix_file_path(self):
        s = "bi_D_contract_split_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, BiColoredHairyGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_contract_split_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, BiColoredHairyGraphComplex.graph_type, self.sub_type, s)


class VertexLoopDegSlice(GraphVectorSpace.DegSlice):
    def __init__(self, deg, h_a_min, h_b_min, even_edges, even_hairs_a, even_hairs_b):
        self.h_a_min = h_a_min
        self.h_b_min = h_b_min
        super(VertexLoopDegSlice, self).__init__(
            [BiColoredHairyGraphComplex.BiColoredHairyGraphVS(v, deg - v, self.h_a_min + v, self.h_b_min + v,
                                                              even_edges, even_hairs_a, even_hairs_b)
             for v in range(0, deg + 1)], deg)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('deg', self.deg), ('min_hairs_a', self.h_a_min), ('min_hairs_b', self.h_b_min)])

    def __eq__(self, other):
        return self.deg == other.deg and self.h_a_min == other.h_a_min and self.h_b_min == other.h_b_min


class VertexLoopBigradedSumVS(GraphVectorSpace.SumVectorSpace):
    def __init__(self, deg_range, h_a_min_range, h_b_min_range, even_edges, even_hairs_a, even_hairs_b):
        self.deg_range = deg_range
        self.h_a_min_range = h_a_min_range
        self.h_b_min_range = h_b_min_range
        self.sub_type = BiColoredHairyGraphComplex.get_sub_type(even_edges, even_hairs_a, even_hairs_b)
        max_deg = max(self.deg_range)
        super(VertexLoopBigradedSumVS, self).__init__(
            [VertexLoopDegSlice(deg, h_a_min + (max_deg - deg), h_b_min + (max_deg - deg), even_edges, even_hairs_a,
                                even_hairs_b) for (deg, h_a_min, h_b_min) in
             itertools.product(self.deg_range, self.h_a_min_range, self.h_b_min_range)])

    def get_type(self):
        return '%s graphs with %s' % (BiColoredHairyGraphComplex.graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('deg', self.deg_range), ('min_hairs_a', self.h_a_min_range),
                                   ('min_hairs_b', self.h_b_min_range)])


class ContractSplitD(GraphOperator.Differential):
    def __init__(self, graded_sum_vs):
        super(ContractSplitD, self).__init__(graded_sum_vs, ContractSplitBiOM.generate_op_matrix_list(graded_sum_vs))

    def get_type(self):
        return 'contract edges and split edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_split_edges_D_%s_%s" % (BiColoredHairyGraphComplex.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, BiColoredHairyGraphComplex.graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        deg_range = self.sum_vector_space.deg_range
        h_a_min_min = min(self.sum_vector_space.h_a_min_range)
        h_a_min_max = max(self.sum_vector_space.h_a_min_range) + (max(deg_range) - min(deg_range))
        h_a_min_range = range(h_a_min_min, h_a_min_max + 1)
        h_b_min_min = min(self.sum_vector_space.h_a_min_range)
        h_b_min_max = max(self.sum_vector_space.h_a_min_range) + (max(deg_range) - min(deg_range))
        h_b_min_range = range(h_b_min_min, h_b_min_max + 1)
        return Shared.OrderedDict([('deg', deg_range), ('min_hairs_a', h_a_min_range),
                                   ('min_hairs_b', h_b_min_range)])

    def get_cohomology_plot_parameter_order(self):
        return (1, 2, 0)


class BiColoredHaryContractSplitBiGC(GraphComplex.GraphComplex):
    def __init__(self, deg_range, h_a_min_range, h_b_min_range, even_edges, even_hairs_a, even_hairs_b):
        self.sub_type = BiColoredHairyGraphComplex.get_sub_type(even_edges, even_hairs_a, even_hairs_b)

        graded_sum_vs = VertexLoopBigradedSumVS(deg_range, h_a_min_range, h_b_min_range, even_edges, even_hairs_a,
                                                even_hairs_b)
        super(BiColoredHaryContractSplitBiGC, self).__init__(graded_sum_vs, [ContractSplitD(graded_sum_vs)])

    def __str__(self):
        return '<%s graphs bi-complex with %s>' % (BiColoredHairyGraphComplex.graph_type, str(self.sub_type))