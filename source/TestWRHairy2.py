import unittest
import itertools
import logging
import Log
import TestGraphComplex
import WRHairyGraphComplex
from sage.all import *


tt = WRHairyGraphComplex.WRHairyGraphVS(4, 3, 1, 1)
ss = WRHairyGraphComplex.WRHairyGraphVS(5, 3, 1, 1)
uu = WRHairyGraphComplex.WRHairyGraphVS(3, 3, 1, 1)

dd1 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(4, 3, 1, 1)
dd2 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(5, 3, 1, 1)

D1 = dd1.get_matrix()
D2 = dd2.get_matrix()
# tt.build_basis()
# tt.plot_all_graphs_to_file(skip_existing=True)
# tt.display_basis_plots()
D1 = D1.dense_matrix().change_ring(QQ)
D2 = D2.dense_matrix().change_ring(QQ)

print(D1)
print(D2)

print("Ranks: ", D1.rank(), D2.rank())
print(ss.get_dimension(), "-D2->", tt.get_dimension(), "-D1->", uu.get_dimension())


nsp = D1.right_kernel()  # (basis='pivot')
print(nsp)
img = D2.transpose().image()
print(img)

HH = nsp/img
print(HH)

nsp2 = D2.kernel()
img2 = D1.image()

print(nsp2)
print(img2)
