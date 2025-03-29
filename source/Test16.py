import OrdinaryGraphComplex 
import SpecialGraphs
import RepresentativeCheck
import GraphOperator
import Parameters
import os


class DualRepresentativeCheck_loopwheel(GraphOperator.OperatorMatrix):
    """Class and method for checking whether a cocycle represents a non-zero cohomology class.
    """

    def __init__(self, k):
        self.op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(8*k+2,4*k+2,False)
        self.V1 = self.op.get_domain()
        self.V2 = self.op.get_target()
        self.k = k
        self.sub_type = OrdinaryGraphComplex.sub_types[False]

        super().__init__(self.V1, RepresentativeCheck.DummyVSOneMore(self.V2))

    def __str__(self):
        return "DualRepresentativeCheck_loopwheel "+str(self.k) 

    def get_matrix_file_path(self):
        s = f"DualRepresentativeCheck_loopwheel{str(self.k)}.txt"
        return os.path.join(Parameters.data_dir, OrdinaryGraphComplex.graph_type, self.sub_type, s)

    def get_rank_file_path(self):
        s = f"DualRepresentativeCheck_loopwheel{str(self.k)}_rank.txt"
        return os.path.join(Parameters.data_dir, OrdinaryGraphComplex.graph_type, self.sub_type, s)

    def is_valid(self):
        return self.op.is_valid()

    @staticmethod
    def is_match(domain, target):
        return True

    def generate_vector(self):
        """returns the vector of the cocycle as a list of pairs (G, coeff)"""
        return SpecialGraphs.doubleloop_graphs(self.k, aut_normalize=True)
        

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        """concatenates the vector to the matrix of op2"""
        if not self.op.is_valid():
            return
        if (not ignore_existing_files) and self.exists_matrix_file():
            return

        vlst = self.V1.graph_list_to_vector(self.generate_vector())
        (mlist, (m,n)) = self.op._load_matrix_list()
        # print("*******", m,n,mlist)
        # newshape = (m+1,n)
        # newmlist = mlist + [(m,k,val) for (k,val) in enumerate(vlst) if val != 0 ]
        newshape = (m,n+1)
        newmlist = mlist + [(k,n,val) for (k,val) in enumerate(vlst) if val != 0 ]


        self._store_matrix_list(newmlist, newshape)





# l=6
# v=10
# k=1

l=10
v=18
k=2

# D = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v,l,False)
# V1 = D.get_domain()
# print("1")
# vv = SpecialGraphs.doubleloop_graphs(k)
# print("2")
# vv2 = SpecialGraphs.doubleloop_graphs(k, aut_normalize=True)
# # A = D.get_matrix()
# print("3")

# v = V1.graph_list_to_vector(vv)
# print("4")

# v2 = V1.graph_list_to_vector(vv2)
# print("5")
# print(v)
# # print(A) 
# print(v2)
k=2
OO = DualRepresentativeCheck_loopwheel(k)
OO.build_matrix(ignore_existing_files=True)
OO.compute_rank(sage="integer", ignore_existing_files=True)
OO.op.compute_rank(sage="integer")
print("New rank: " ,OO.get_matrix_rank())
# old rank
print("Dimension: ", OO.get_domain().get_dimension())
print("Old Rank: ", OO.op.get_matrix_rank())