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


class ForestedMoritaTetrahedronCheck(RepresentativeCheck.RepresentativeCheck):
    # TODO: not finished yet...
    # The forested wheel with 2k+1 spokes and all but one of the edges along the rim marked
    def __init__(self, k):
        even_edges = True
        self.even_edges = even_edges
        self.k = k
        self.sub_type = ForestedGraphComplex.sub_types[self.even_edges]
        op1 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(2*k+1, 2*k+2, 1, even_edges)
        op2 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(2*k+1, 2*k+1, 1, even_edges)
        name = f"Forested wheel class "+self.sub_type
        super(ForestedMoritaTetrahedronCheck, self).__init__(op1, op2, name)

    def get_matrix_file_path(self):
        s = f"contract_unmarkD_forested_wheel_check.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = f"contract_unmarkD_forested_wheel_check_rank.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def generate_vector(self):
        k = self.k
        G = Graph(4*k + 6) # 2k+3 normal vertices then 2k+2 bivalent edge vertices, then 1 hair
        # the first vertex is the center, then the spokes-rim-vertices, then the special rim vertex
        # attach spokes 
        for j in range(2*k+1):
            G.add_edge(0, 2*k+3+j)
            G.add_edge(j+1, 2*k+3+j)
        # edge from hair to special rim vertex
        G.add_edge(2*k+4, 4*k+5)

        ret = []

        # attach all 2k+2 rim edges... all but one is marked
        for j in range(2*k+2):
            # j is the edge that is not marked
            #
            pass

        return ret 

class ForestedWheelCheck(RepresentativeCheck.RepresentativeCheck):
    def __init__(self, even_edges):
        self.even_edges = even_edges
        self.sub_type = ForestedGraphComplex.sub_types[even_edges]
        op1 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7, 8, 0, even_edges)
        op2 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7, 8, 0, even_edges)
        name = f"Forested Morita tetrahedron class "+self.sub_type
        super(ForestedMoritaTetrahedronCheck, self).__init__(op1, op2, name)

    def get_matrix_file_path(self):
        s = f"contract_unmarkD_morita_tetrahedron_check.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = f"contract_unmarkD_morita_tetrahedron_check_rank.txt"
        return os.path.join(Parameters.data_dir, ForestedGraphComplex.graph_type, self.sub_type, s)

    def generate_vector(self):
        return [ (SpecialGraphs.forested_morita_tetrahedron([p1,p2,p3,p4]),
                     Shared.Perm(p1).signature() * Shared.Perm(p2).signature() 
                     * Shared.Perm(p3).signature() * Shared.Perm(p4).signature()) 
                    for p1 in itertools.permutations(range(3))
                    for p2 in itertools.permutations(range(3))
                    for p3 in itertools.permutations(range(3))
                    for p4 in itertools.permutations(range(3))
                    ]

##### run tests #####
# OrdinaryWheelCheck(3, False).checkit(sage="integer")
# OrdinaryWheelCheck(4, False).checkit(sage="integer")
# OrdinaryWheelCheck(5, False).checkit(sage="integer")

# OrdinaryWheelCheck(3, True).checkit(sage="integer")
# OrdinaryWheelCheck(4, True).checkit(sage="integer")
# OrdinaryWheelCheck(5, True).checkit(sage="integer")


# ForestedRingCheck(3, True).checkit(linbox="mod")
# ForestedRingCheck(4, True).checkit(linbox="mod")
# ForestedRingCheck(5, True).checkit(linbox="mod")

# ForestedMoritaCheck(3, False).checkit(linbox="mod")
# ForestedMoritaCheck(4, False).checkit(linbox="mod")
# ForestedMoritaCheck(5, False).checkit(linbox="mod")


FC = ForestedMoritaTetrahedronCheck(False)
# vvv=FC.generate_vector()
# (G,a) = vvv[0]
# G.show()
# print(G.graph6_string())
# (G,a) = vvv[3]
# print(G.graph6_string())
# G.show()

# print(FC.is_cocycle())
print("matrix...")
# FC.build_matrix()
print("built.. starting rank")
# FC.compute_rank(linbox="modprecond")
print("Done.")
FC.is_cocycle_exact()