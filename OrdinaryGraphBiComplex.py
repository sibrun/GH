from sage.all import *
import itertools
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH
import Parameters
import OrdinaryGraphComplex as OGC


class CeDeleBiOM(GO.BiOperatorMatrix):
    def __init__(self, domain, target):
        self.sub_type = domain.get_vs_list()[0].sub_type
        super(CeDeleBiOM, self).__init__(domain, target, OGC.ContractEdgesGO, OGC.DeleteEdgesGO)

    @staticmethod
    def is_match(domain, target):
        (d_deg,) = domain.get_ordered_param_dict().get_value_tuple()
        (t_deg,) = target.get_ordered_param_dict().get_value_tuple()
        return d_deg - 1 == t_deg

    def get_matrix_file_path(self):
        s = "bi_D_ce_dele_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, OGC.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_ce_dele_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, OGC.graph_type, self.sub_type, s)

    def get_ref_matrix_file_path(self):
        pass

    def get_ref_rank_file_path(self):
        pass


class VertexLoopDegSlice(GVS.DegSlice):
    def __init__(self, deg, even_edges):
        self.even_edges = even_edges
        super(VertexLoopDegSlice, self).__init__(
            [OGC.OrdinaryGVS(v, deg - v, self.even_edges) for v in range(0, deg + 1)], deg)

    def get_ordered_param_dict(self):
        return SH.OrderedDict([('deg', self.deg)])

    def __eq__(self, other):
        return self.deg == other.deg


class VertexLoopBigradedSumVS(GVS.SumVectorSpace):
    def __init__(self, deg_range, even_edges):
        self.deg_range = deg_range
        self.even_edges = even_edges
        self.sub_type = OGC.sub_types.get(even_edges)
        super(VertexLoopBigradedSumVS, self).__init__([VertexLoopDegSlice(deg, even_edges) for deg in deg_range])

    def get_type(self):
        return '%s graphs with %s' % (OGC.graph_type, self.sub_type)

    def get_ordered_param_range_dict(self):
        return SH.OrderedDict([('deg', self.deg_range)])


class CeDeleD(GO.Differential):
    def __init__(self, graded_sum_vs):
        super(CeDeleD, self).__init__(graded_sum_vs, CeDeleBiOM.generate_op_matrix_list(graded_sum_vs))

    def get_type(self):
        return 'contract edges and delete edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_edges_delete_edges_D_%s_%s.png" % (OGC.graph_type, sub_type)
        return os.path.join(Parameters.plots_dir, OGC.graph_type, sub_type, s)

    def get_ordered_cohomology_param_range_dict(self):
        deg_range = self.sum_vector_space.deg_range
        return SH.OrderedDict([('deg', deg_range)])

    def get_cohomology_plot_parameter_order(self):
        return (0, 1)


class OrdinaryCeDeleBiGC(GC.GraphComplex):
    def __init__(self, deg_range, even_edges):
        self.deg_range = deg_range
        self.even_edges = even_edges
        self.sub_type = OGC.sub_types.get(self.even_edges)
        graded_sum_vs = VertexLoopBigradedSumVS(deg_range, even_edges)
        super(OrdinaryCeDeleBiGC, self).__init__(graded_sum_vs, [CeDeleD(graded_sum_vs)])

    def __str__(self):
        return '<ordinary graphs bi-complex with %s>' % str(self.sub_type)