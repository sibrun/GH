import unittest
import itertools
import logging
import Log
import TestGraphComplex
import OrdinaryGraphComplex
from sage.all import *
import StoreLoad
import OrdinaryMerkulovComplex
import GraphOperator
import ForestedGraphComplex
import os


# FD = ForestedGraphComplex.ForestedTopDegSlice(4, 3, 1, True, 1)
# FD.build_basis()
# FD.display_basis_plots()

# FGC = ForestedGraphComplex.ForestedGC(range(15),
#                                       range(5), range(6), range(3), True, {'contract'})
# FGC.build_basis()


op1 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
    4, 3, 2, True)
vs = op1.get_domain()
# vs.display_basis_plots()
# op2 = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
#     4, 4, 1, True)
op1.build_matrix()
# op2.build_matrix()

D = op1.get_matrix()
# DD = op2.get_matrix()

# op3 = ForestedGraphComplex.ContractEdgesGO.generate_operator(5, 3, 3, 1, True)
# op3.build_matrix()
# Dc = op3.get_matrix()

print(D)

# print(DD)
# print(Dc)

kk = D.right_kernel()
print(kk)

# hunt smallest element
minlen = 100000
# min(len(w.nonzero_positions()) for w in kk.basis())
for v in kk.basis():
    nz = v.nonzero_positions()
    if len(nz) <= minlen+2:
        print(nz)
        if len(nz) < minlen:
            minlen = len(nz)
