import unittest
import itertools
import logging
import Log
import TestGraphComplex
import ForestedGraphComplex
from sage.all import *

# VS = ForestedGraphComplex.ForestedGVS(6, 4, 3, 0, True)
# VS = ForestedGraphComplex.ForestedDegSlice(4, 3, 0, True)
# VS.build_basis()
# VS.display_basis_plots()

# FTD = ForestedGraphComplex.ContractUnmarkTopD(
#     range(6), range(15), range(1), True)

# FTD.build_matrix()
# FTD.compute_rank(sage="integer")

# print(FTD.get_cohomology_dim_dict())

# FTD.plot_cohomology_dim()


FD = ForestedGraphComplex.ForestedDegSlice(2,0,2, True)
print(FD.is_valid())
pop = FD.get_isotypical_projector(0)
pop.build_matrix()
