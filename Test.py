import OrdinaryGraphComplex as OGC
from sage.all import *
import Parameters

if __name__ == "__main__":

    op = OGC.ContractEdgesGO.generate_operator(10,9,True)
    M = op.get_matrix_transposed(32633)

