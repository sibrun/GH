from sage.all import *
import OrdinaryGraphComplex
import random


op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(13,8,False)
D = op.get_matrix()
if D.nrows() > D.ncols():
    D = D.transpose()
r = D.rank()
# select r rows at random
rows = random.sample(range(D.nrows()), min(r+0, D.nrows()))
D_sub = D[rows, :]
r_sub = D_sub.rank()
print("Size of D:", D.nrows(), "x", D.ncols())
print("Original rank:", r)
print("Submatrix rank:", r_sub)