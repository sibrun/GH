import ForestedGraphComplex
import CHairyGraphComplex
from sage.all import *

from timeit import default_timer as timer

vs0 = ForestedGraphComplex.ForestedDegSlice(5,6,1,True)
vs1 = ForestedGraphComplex.ForestedDegSlice(5,5,1,True)
vs2 = ForestedGraphComplex.ForestedDegSlice(5,4,1,True)

# vs1.display_basis_plots()

op1 = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(5,6,1,True)
op2 = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(5,5,1,True)

vs0.build_basis()
vs1.build_basis()
vs2.build_basis()
op1.build_matrix()
op2.build_matrix()
op1.compute_rank(sage="integer")
op2.compute_rank(sage="integer")

print( vs0.get_dimension(),vs1.get_dimension(),vs2.get_dimension() )

print(op1.get_matrix_rank(), op2.get_matrix_rank())

D1 = op1.get_matrix()
D2 = op2.get_matrix()
print(D1)
print(D2)

# v0=Matrix(24,1)
# v0[15] = 1
# v0[18] = -1
# print(v0)
# print(D2 * v0)

# v1=Matrix(24,1)
# v1[5] = 1
# v1[16] = 1
# print(v1)
# print(D2 * v1)
# print(D1.rank(), D1.augment(v0).rank(), D1.augment(v1).rank(), D1.augment(v0).augment(v1).rank())

# print(op2.get_matrix().right_kernel())
# print(op1.get_matrix().transpose().image())

