from sage.all import *
import itertools
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH
import Parameters
import HairyGraphComplex as HGC


class CeEt1hBiOM(GO.BiOperatorMatrix):
    def __init__(self, domain, target, operator_cls1, operator_cls2):
        self.sub_type = domain.get_vs_list()[0].sub_type
        super(CeEt1hBiOM, self).__init__(domain, target, operator_cls1, operator_cls2)

    @staticmethod
    def is_match(domain, target):
        (d_deg, d_h_min) = domain.get_ordered_param_dict().get_value_tuple()
        (t_deg, t_h_min) = target.get_ordered_param_dict().get_value_tuple()
        return (d_deg, d_h_min) == (t_deg + 1, t_h_min - 1)

    def get_matrix_file_path(self):
        s = "bi_D_ce_et1h_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, HGC.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_ce_et1h_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, HGC.graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        pass

    def get_ref_rank_file_path(self):
        pass


class VertexLoopDegSlice(GVS.DegSlice):
    def __init__(self, deg, h_min, even_edges, even_hairs):
        self.h_min = h_min
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        super(VertexLoopDegSlice, self).__init__(
            [HGC.HairyGraphVS(v, deg - v, self.h_min + v, self.even_edges, self.even_hairs)
             for v in range(0, deg + 1)], deg)

    def get_ordered_param_dict(self):
        return SH.OrderedDict([('deg', self.deg), ('min_hairs', self.h_min)])

    def __eq__(self, other):
        return self.deg == other.deg and self.h_min == other.h_min


class VertexLoopBigradedSumVS(GVS.SumVectorSpace):
    def __init__(self, deg_range, h_min_range, even_edges, even_hairs):
        self.deg_range = deg_range
        self.h_min_range = h_min_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = HGC.sub_types.get((even_edges, even_hairs))
        max_deg = max(self.deg_range)
        super(VertexLoopBigradedSumVS, self).__init__(
            [VertexLoopDegSlice(deg, h_min + (max_deg - deg), even_edges, even_hairs) for (deg, h_min) in
             itertools.product(self.deg_range, self.h_min_range)])

    def get_type(self):
        return '%s graphs with %s' % (HGC.graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return SH.OrderedDict([('deg', self.deg_range), ('min_hairs', self.h_min_range)])


class CeEt1hD(GO.BiDifferential):
    def __init__(self, graded_sum_vs):
        super(CeEt1hD, self).__init__(graded_sum_vs, HGC.ContractEdgesGO, HGC.EdgeToOneHairGO, CeEt1hBiOM)

    def get_type(self):
        return 'contract edges and edge to one hair'

    def get_cohomology_plot_path(self):
        sub_type = self.graded_sum_vs.sub_type
        s = "cohomology_dim_contract_edges_edge_to_one_hair_D_%s_%s.png" % (HGC.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, HGC.graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        deg_range = self.graded_sum_vs.deg_range
        h_min_min = min(self.graded_sum_vs.h_min_range)
        h_min_max = max(self.graded_sum_vs.h_min_range) + (max(deg_range) - min(deg_range))
        h_range = range(h_min_min, h_min_max + 1)
        return SH.OrderedDict([('deg', deg_range), ('min_hairs', h_range)])

    def get_cohomology_plot_parameter_order(self):
        return (0, 1)


class HairyBiGC(GC.GraphComplex):
    def __init__(self, deg_range, h_min_range, even_edges, even_hairs):
        self.deg_range = deg_range
        self.h_min_range = h_min_range
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        graded_sum_vs = VertexLoopBigradedSumVS(deg_range, h_min_range, even_edges, even_hairs)
        super(HairyBiGC, self).__init__(graded_sum_vs, [CeEt1hD(graded_sum_vs)])

    def __str__(self):
        return '<hairy graphs bi-complex with %s>' % str(self.sub_type)

