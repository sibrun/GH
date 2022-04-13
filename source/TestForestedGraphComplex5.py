import unittest
import itertools
import logging
import Log
import TestGraphComplex
import ForestedGraphComplex
from sage.all import *


max_loops = 4
maxh = 0
max_verts = 2*max_loops - 2 + maxh
maxm = max_verts+1

# PFGC = ForestedGraphComplex.PreForestedGraphSumVS(
#     range(0, max_verts+1), range(0, max_loops+1), range(0, maxm+1), range(0, max_loops+1+maxh))
# PFGC.build_basis(ignore_existing_files=False)


evenedges = True

# FGC = ForestedGraphComplex.ForestedGC(
#     range(0, max_verts+1), range(0, max_loops+1), range(0, maxm+1), range(0, maxh+1), evenedges, {'contract', 'unmark'})
# FGC.build_basis(ignore_existing_files=False)
# FGC.build_matrix(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=False)

# FBGC = ForestedGraphComplex.ForestedContractUnmarkBiGC(
#     range(0, max_loops+1), range(0, maxm+1), range(0, maxh+1), evenedges)
# FBGC.build_basis(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=False)
# FBGC.build_matrix(progress_bar=False, info_tracker=False,
#                   ignore_existing_files=False)
# # FBGC.square_zero_test()

# FBGC.compute_rank(ignore_existing_files=False, sage="integer")


# FBGC.print_dim_and_eulerchar()
# FBGC.print_cohomology_dim()

Op = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
    4, 3, 0, evenedges)
vs1 = Op.get_domain()
vs2 = Op.get_target()

# vs1.display_basis_plots()

M = Op.get_matrix()

print(M)
