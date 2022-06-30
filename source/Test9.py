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


vs = ForestedGraphComplex.ForestedGVS(2,1,1,2,True)

opP0 = vs.get_isotypical_projector(0)
opP1 = vs.get_isotypical_projector(1)

opP0.build_matrix(ignore_existing_files=True)
opP1.build_matrix(ignore_existing_files=True)

# for opP in [opP1]:
for opP in [opP0, opP1]:
    print(opP, opP.rep_index)
    b = vs.get_basis()
    for G in b:
        print(vs.graph_to_canon_g6(G))
        for (GG, v) in opP.operate_on(G):
            print(vs.graph_to_canon_g6(GG), v)
        

    