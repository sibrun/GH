import unittest
import itertools
import logging
import Log
import TestGraphComplex
import ForestedGraphComplex
from sage.all import *


maxl = 4
maxh = 0
maxv = 2*maxl - 2 + maxh
maxm = maxv+1

PFGC = ForestedGraphComplex.PreForestedGraphSumVS(range(0,maxv+1), range(0,maxl+1), range(0,maxm+1), range(0,maxh+1))
PFGC.build_basis(ignore_existing_files=True)


evenedges = True

FGC = ForestedGraphComplex.ForestedGC(
    range(0, maxv+1), range(0, maxl+1), range(0, maxm+1), range(0, maxh+1), evenedges, {'contract', 'unmark'})
FGC.build_basis(ignore_existing_files=True)
FGC.build_matrix(progress_bar=False, info_tracker=False,
                 ignore_existing_files=True)
# FGC.square_zero_test()

FBGC = ForestedGraphComplex.ForestedContractUnmarkBiGC(
    range(0, maxl+1), range(0, maxm+1), range(0, maxh+1), evenedges)
FBGC.build_basis(progress_bar=False, info_tracker=False,
                 ignore_existing_files=True)
FBGC.build_matrix(progress_bar=False, info_tracker=False,
                  ignore_existing_files=True)
FBGC.square_zero_test()

FBGC.compute_rank(ignore_existing_files=False, sage="integer")



FBGC.print_dim_and_eulerchar()
FBGC.print_cohomology_dim()

uc1 = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(4,3,0,evenedges)
uc2 = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(4,4,0,evenedges)

D1 = uc1.get_matrix()
D2 = uc2.get_matrix()
print(D1*D2) # must be=0

print(D1) # for info
print(D2)


uc1.domain.display_basis_plots()
g_ind = 42 #index of graph in basis
v=vector( [0 for j in range(43)] )
v[g_ind]=1
print(v)

print(D1*v) # should be 0

