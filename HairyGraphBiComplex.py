import itertools
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH
import NautyInterface as NI
import OrdinaryGraphComplex as OGC
import StoreLoad as SL
import Parameters
import HairyGraphComplex as HGC


v_range = range(0,8)
l_range = range(0,8)
h_range = range(0,10)
even_edges = False
even_hairs = False

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


class CeEt1hD(GO.Differential):
    def __init__(self, graded_vs, dif_1, dif_2):
        super(CeEt1hD, self).__init__(vector_space, CeEt1hBiOM.generate_op_matrix_list(graded_vs, dif_1, dif_2))

    def get_type(self):
        return 'contract edges and edge to one hair'


class EdgeToOneHairGC(GC.GraphComplex):
    def __init__(self, vector_space, differential):
        self.sub_type = HGC.sub_types.get(even_edges, even_hairs)
        super(EdgeToOneHairGC, self).__init__(vector_space, differential)

    def get_cohomology_plot_path(self):
        s = "cohomology_dim_edge_to_one_hair_%s_%s.png" % (HGC.graph_type, self.sub_type)
        return os.path.join(Parameters.plots_dir, HGC.graph_type, self.sub_type, s)


def generate_deg_slice(deg, h_min):
    return GVS.DegSlice([HGC.HairyGVS(n, deg - n, h_min) for n in range(deg)], deg)


vector_space = HGC.SumHairyGVS(v_range, l_range, h_range, even_edges, even_hairs)
dif_1 = HGC.ContractEdgesD(vector_space)
dif_2 = HGC.EdgeToOneHairD(vector_space)

graded_vs = GVS.SumVectorSpace([generate_deg_slice(4 + n, n) for n in range(0, 3)])

dif = CeEt1hD(graded_vs, dif_1, dif_2)






