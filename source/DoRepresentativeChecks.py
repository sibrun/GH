import RepresentativeCheck
import OrdinaryGraphComplex
import ForestedGraphComplex
import Parameters
import SpecialGraphs
import os
import itertools
import Shared

class OrdinaryWheelCheck(RepresentativeCheck.RepresentativeCheck):
    def __init__(self, k_spokes, even_edges):
        self.k_spokes = k_spokes 
        self.even_edges = even_edges
        self.sub_type = OrdinaryGraphComplex.sub_types[even_edges]
        op1 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(k_spokes+1, k_spokes, even_edges)
        op2 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(k_spokes+2, k_spokes, even_edges)
        name = f"Wheel ordinary ({k_spokes} spokes) "+self.sub_type
        super(OrdinaryWheelCheck, self).__init__(op1, op2, name)

    def get_matrix_file_path(self):
        s = f"contractD_wheelcheck{self.k_spokes}.txt"
        return os.path.join(Parameters.data_dir, OrdinaryGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = f"contractD_wheelcheck{self.k_spokes}_rank.txt"
        return os.path.join(Parameters.data_dir, OrdinaryGraphComplex.graph_type, self.sub_type, s)

    def generate_vector(self):
        return [ (SpecialGraphs.wheel_graph(self.k_spokes), 1) ]


class ForestedRingCheck(RepresentativeCheck.RepresentativeCheck):
    def __init__(self, n_marked, even_edges):
        self.n_marked = n_marked 
        self.even_edges = even_edges
        self.sub_type = ForestedGraphComplex.sub_types[even_edges]
        op1 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(n_marked+1, n_marked, 0, even_edges)
        op2 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(n_marked+1, n_marked+1, 0, even_edges)
        name = f"Forested ring ({n_marked} marked edges) "+self.sub_type
        super(ForestedRingCheck, self).__init__(op1, op2, name)

    def get_matrix_file_path(self):
        s = f"contract_unmarkD_ringcheck{self.n_marked}.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = f"contract_unmarkD_ringcheck{self.n_marked}_rank.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def generate_vector(self):
        return [ (SpecialGraphs.forested_ring_graph(self.n_marked), 1) ]


class ForestedMoritaCheck(RepresentativeCheck.RepresentativeCheck):
    def __init__(self, k, even_edges):
        self.k = k
        self.even_edges = even_edges
        self.sub_type = ForestedGraphComplex.sub_types[even_edges]
        op1 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(k+1, 2*k-2, 0, even_edges)
        op2 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(k+1, 2*k-1, 0, even_edges)
        name = f"Forested Morita class ({k} bridge edges) "+self.sub_type
        super(ForestedMoritaCheck, self).__init__(op1, op2, name)

    def get_matrix_file_path(self):
        s = f"contract_unmarkD_moritacheck{self.k}.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = f"contract_unmarkD_moritacheck{self.k}_rank.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def generate_vector(self):
        return [ (SpecialGraphs.forested_morita_graph(self.k, p), Shared.Perm(p).signature()) 
                    for p in itertools.permutations(range(self.k))]



##### run tests #####
# OrdinaryWheelCheck(3, False).checkit(sage="integer")
# OrdinaryWheelCheck(4, False).checkit(sage="integer")
# OrdinaryWheelCheck(5, False).checkit(sage="integer")

# OrdinaryWheelCheck(3, True).checkit(sage="integer")
# OrdinaryWheelCheck(4, True).checkit(sage="integer")
# OrdinaryWheelCheck(5, True).checkit(sage="integer")


ForestedRingCheck(3, True).checkit(linbox="mod")
ForestedRingCheck(4, True).checkit(linbox="mod")
ForestedRingCheck(5, True).checkit(linbox="mod")

ForestedMoritaCheck(3, False).checkit(linbox="mod")
ForestedMoritaCheck(4, False).checkit(linbox="mod")
ForestedMoritaCheck(5, False).checkit(linbox="mod")