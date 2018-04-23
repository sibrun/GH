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



def generate_deg_slice(deg, h_min):
    return GVS.DegSlice([HGC.HairyGVS(n, deg - n, h_min) for n in range(deg)], deg)

GradedVS = GVS.SumVectorSpace([generate_deg_slice(5 + n, n) for n in range(0, 3)])





