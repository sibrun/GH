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

    diff = D1*D2
    diffs = sum(abs(c) for cc in diff.columns() for c in cc)
    print(diffs)

    isocomp_dim = P1.rank()
    r1 = D1P.rank()
    r2 = D2P.rank()
    print(isocomp_dim, r1, r2)
    cohomdim = isocomp_dim - r1-r2
    part = Partitions(n_hairs)[rep_ind]
    rep_dim = symmetrica.charvalue(part, [1 for j in range(n_hairs)])
    if cohomdim > 0:
        print("Cohomology found:  w=", n_ws, ", h=", n_hairs, ", l=", n_loops, ", vertices=", n_vertices,
              " (degree ", n_vertices+n_loops +
              1, "), partition=", part,  ", invpartition=", part.conjugate(),
              ", multiplicity=", cohomdim/rep_dim, ", cohomdim=", cohomdim)
    return isocomp_dim - r1-r2


def getCohomDimPAll(gvs):
    for w in gvs.w_range:
        for h in gvs.h_range:
            for l in gvs.l_range:
                for v in gvs.v_range:
                    D1 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
                        v, l, h, w)
                    D2 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
                        v+1, l, h, w)
                    try:
                        d = WRHairyGraphComplex.WRHairyGraphVS(
                            v, l, h, w).get_dimension()
                        r1 = D1.get_matrix_rank()
                        r2 = D2.get_matrix_rank()
                        if d-r1-r2 > 0:
                            for rep_ind in range(len(Partitions(h))):
                                getCohomDimP(v, l, h, w, rep_ind)
                    except:
                        pass


WGC = WRHairyGraphComplex.WRHairyGC(range(14), range(
    3, 4), range(4, 5), range(1, 2), ['contract'])

WGC.build_basis(progress_bar=False, info_tracker=False,
                ignore_existing_files=True)
WGC.build_matrix(progress_bar=False, info_tracker=False,
                 ignore_existing_files=True)

# # WGC.build_basis(progress_bar=False, info_tracker=False, ignore_existing_files=False)
# # WGC.build_matrix(progress_bar=False, info_tracker=False, ignore_existing_files=False)

# # WGC.square_zero_test()

# # WGC.compute_rank(ignore_existing_files=True, sage="mod")
WGC.compute_rank(ignore_existing_files=True, sage="integer")
# # WGC.plot_cohomology_dim(to_html=True)
# # Euler char
WGC.print_dim_and_eulerchar()
WGC.print_cohomology_dim()

getCohomDimPAll(WGC)

# print(getCohomDimP(9, 4, 3, 1, 0))
# print(getCohomDimP(9, 4, 3, 1, 1))
# print(getCohomDimP(9, 4, 3, 1, 2))
