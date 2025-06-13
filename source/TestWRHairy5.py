import unittest
import itertools
import logging
import Log
import TestGraphComplex
import WRHairyGraphComplex
from sage.all import *


def getCohomDimP(n_vertices, n_loops, n_hairs, n_ws, rep_ind):
    tt = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
        n_vertices, n_loops, n_hairs, n_ws)
    tu = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
        n_vertices+1, n_loops, n_hairs, n_ws)
    symmp1 = WRHairyGraphComplex.SymmProjector.generate_operator(
        n_vertices, n_loops, n_hairs, n_ws, rep_ind)

    D1 = tt.get_matrix()
    D2 = tu.get_matrix()
    # C = D2*D1
    symmp1.build_matrix(ignore_existing_files=True)
    P1 = symmp1.get_matrix()
    print("matrices loaded")

    D1P = D1*P1
    D2P = P1*D2
    print("computing ranks....")

    # diff = D1*D2
    # diffs = sum(abs(c) for cc in diff.columns() for c in cc)
    # print(diffs)

    isocomp_dim = P1.rank()
    r1 = D1P.rank()
    r2 = D2P.rank()
    print(isocomp_dim, r1, r2)
    return isocomp_dim - r1-r2


WGC = WRHairyGraphComplex.WRHairyGC(range(14), range(
    3, 4), range(4, 5), range(2, 3), ['contract'])

# WGC.build_basis(progress_bar=False, info_tracker=False,
#                 ignore_existing_files=True)
# WGC.build_matrix(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=True)

# # # WGC.build_basis(progress_bar=False, info_tracker=False, ignore_existing_files=False)
# # # WGC.build_matrix(progress_bar=False, info_tracker=False, ignore_existing_files=False)

# # # WGC.square_zero_test()

# # # WGC.compute_rank(ignore_existing_files=True, sage="mod")
# WGC.compute_rank(ignore_existing_files=True, sage="integer")
# # # WGC.plot_cohomology_dim(to_html=True)
# # # Euler char
# WGC.print_dim_and_eulerchar()
# WGC.print_cohomology_dim()

print(getCohomDimP(5, 3, 4, 2, 0))
print(getCohomDimP(5, 3, 4, 2, 1))
print(getCohomDimP(5, 3, 4, 2, 2))
print(getCohomDimP(5, 3, 4, 2, 3))
print(getCohomDimP(5, 3, 4, 2, 4))

print(getCohomDimP(6, 3, 4, 2, 0))
print(getCohomDimP(6, 3, 4, 2, 1))
print(getCohomDimP(6, 3, 4, 2, 2))
print(getCohomDimP(6, 3, 4, 2, 3))
print(getCohomDimP(6, 3, 4, 2, 4))

print(getCohomDimP(7, 3, 4, 2, 0))
print(getCohomDimP(7, 3, 4, 2, 1))
print(getCohomDimP(7, 3, 4, 2, 2))
print(getCohomDimP(7, 3, 4, 2, 3))
print(getCohomDimP(7, 3, 4, 2, 4))
