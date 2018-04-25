from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH
import Parameters
import HairyGraphComplex as HGC


class CeEt1hBiOM(GO.BiOperatorMatrix):
    def __init__(self, domain, target, op_collection1, op_collection2):
        self.sub_type = domain.get_vs_list()[0].sub_type
        super(CeEt1hBiOM, self).__init__(domain, target, op_collection1, op_collection2)

    @staticmethod
    def is_match(domain, target):
        return domain.get_deg() == target.get_deg() + 1

    def get_matrix_file_path(self):
        s = "bi_D_ce_et1h_%d.txt" % self.domain.get_deg()
        return os.path.join(Parameters.data_dir, HGC.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = "bi_D_ce_et1h_%d_rank.txt" % self.domain.get_deg()
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
        super(VertexLoopDegSlice, self).__init__([HGC.HairyGVS(n, deg - n, self.h_min + n, self.even_edges, self.even_hairs)
                                                  for n in range(0, deg + 1)], deg)

    def get_ordered_param_dict(self):
        return SH.OrderedDict({'deg': self.deg, 'min_hairs': self.h_min})


class VertexLoopBigradedVS(GVS.SumVectorSpace):
    def __init__(self, deg_range, h_min, even_edges, even_hairs):
        self.deg_range = deg_range
        self.h_min = h_min
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        min_deg = min(self.deg_range)
        super(VertexLoopBigradedVS, self).__init__(
            [VertexLoopDegSlice(n, self.h_min + n - min_deg, even_edges, even_hairs) for n in self.deg_range])

    def get_ordered_param_range_dict(self):
        return SH.OrderedDict({'deg': self.deg_range, 'min_hairs': self.h_min})


class CeEt1hD(GO.BiDifferential):
    def __init__(self, graded_vs):
        super(CeEt1hD, self).__init__(graded_vs, HGC.ContractEdgesD, HGC.EdgeToOneHairD, CeEt1hBiOM)

    def get_type(self):
        return 'contract edges and edge to one hair'


class CeEt1hGC(GC.GraphComplex):
    def __init__(self, deg_range, h_min, even_edges, even_hairs):
        self.deg_range = deg_range
        self.h_min = h_min
        self.even_edges = even_edges
        self.even_hairs = even_hairs
        self.sub_type = HGC.sub_types.get(even_edges, even_hairs)
        graded_vs = VertexLoopBigradedVS(deg_range, h_min, even_edges, even_hairs)
        super(CeEt1hGC, self).__init__(graded_vs,
                                       CeEt1hD(graded_vs))

    def get_cohomology_plot_path(self):
        s = "cohomology_dim_contract_edges_edge_to_one_hair_%s_%s.png" % (HGC.graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, HGC.graph_type, self.sub_type, s)


gc = CeEt1hGC(range(7,9), 0, False, False)
gc.build_basis()
gc.build_matrix()

