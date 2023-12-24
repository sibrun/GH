import unittest
import itertools
import logging
import Log
import TestGraphComplex
import ForestedGraphComplex
from sage.all import *


log_file = "Forested_Unittest.log"


# def check_graphs_vs_basis(GVS, w):
#     # Takes a list of graphs and checks whether they are found in the basis
#     ba = GVS.get_basis_g6()
#     for HH in w:
#         H = Graph(HH) if type(HH) is str else H
#         g6, sgn = GVS.graph_to_canon_g6(H)
#         autom_list = H.automorphism_group(partition=GVS.get_partition()).gens()
#         if GVS._has_odd_automorphisms(H, autom_list):
#             print(g6, " has odd automorphisms")
#         else:
#             if not g6 in ba:
#                 print(g6, " not found in basis")
#             else:
#                 print(g6, " exists with index ", ba.index(g6))


def DSquareTestSingleUnmark(n_vertices, n_loops, n_marked, n_hairs, even_edges, j_to_pick=-1, plot_basis=False):
    tt = ForestedGraphComplex.UnmarkEdgesGO.generate_operator(
        n_vertices, n_loops, n_marked, n_hairs, even_edges)
    tu = ForestedGraphComplex.UnmarkEdgesGO.generate_operator(
        n_vertices, n_loops, n_marked-1, n_hairs, even_edges)
    D1 = tt.get_matrix()
    D2 = tu.get_matrix()
    C = D2*D1

    print(D1)
    print(D2)
    print(C)

    ba0 = tt.domain.get_basis_g6()
    ba1 = tu.domain.get_basis_g6()
    ba2 = tu.target.get_basis_g6()

    if plot_basis:
        tt.domain.display_basis_plots()
        tu.domain.display_basis_plots()
        tu.target.display_basis_plots()

    if (j_to_pick < 0):
        for i in range(C.nrows()):
            for j in range(C.ncols()):
                if C[i, j] != 0:
                    print(i, j, C[i, j])
                    j_to_pick = j
        if j_to_pick < 0:
            print("success, squares to zero")
            return
        else:
            print("Does not square to zero, checking index ",
                  j_to_pick, " g6code ", ba0[j_to_pick])

    G = Graph(ba0[j_to_pick])
    w = tt.operate_on(G)

    # check whether graphs are in basis
    for H, x in w:
        g6, sgn = tu.domain.graph_to_canon_g6(H)
        autom_list = H.automorphism_group(
            partition=tu.domain.get_partition()).gens()
        if tu.domain._has_odd_automorphisms(H, autom_list):
            print(g6, " has odd automorphisms")
        else:
            if g6 not in ba1:
                print(g6, " not found in basis ", " v=", x)
            else:
                print(g6, " exists at index ", ba1.index(g6), " v=", x)

    # compute D^2
    ww = [(HH, x*xx) for H, x in w for HH, xx in tu.operate_on(H)]
    wwd = {}
    for H, x in ww:
        g6, sgn = tu.target.graph_to_canon_g6(H)
        if g6 in wwd:
            wwd[g6] += x
        else:
            wwd[g6] = x
    print(wwd)
    nonzeroflag = false
    for g6, x in wwd.items():
        if x != 0:
            print("Nonzero entry: ", g6, x)
            nonzeroflag = true
    if not nonzeroflag:
        print("all entries zero, i.e., success.")


# def SymmRepDimension(p):
#     return symmetrica.charvalue(p, [1 for j in range(sum(p))])


# def PSquareTest(n_vertices, n_loops, n_hairs, n_ws, rep_ind):
#     # Tests whether the projector squares to itself
#     symmp = WRHairyGraphComplex.SymmProjector.generate_operator(
#         n_vertices, n_loops, n_hairs, n_ws, rep_ind)
#     symmp.build_matrix(ignore_existing_files=False)
#     P = symmp.get_matrix()
#     # should be zero
#     diff = P*P - factorial(n_hairs) * P / \
#         SymmRepDimension(Partitions(n_hairs)[rep_ind])
#     diffs = sum(abs(c) for cc in diff.columns() for c in cc)
#     print("PSquareTest", n_vertices, n_loops,
#           n_hairs, n_ws, rep_ind, ", result: ", diffs)
#     if diffs > 0:
#         print(P)


# def PPTest(n_vertices, n_loops, n_hairs, n_ws, rep_ind1, rep_ind2):
#     # tests whether two projectors have zero product
#     symmp1 = WRHairyGraphComplex.SymmProjector.generate_operator(
#         n_vertices, n_loops, n_hairs, n_ws, rep_ind1)
#     symmp2 = WRHairyGraphComplex.SymmProjector.generate_operator(
#         n_vertices, n_loops, n_hairs, n_ws, rep_ind2)
#     symmp1.build_matrix(ignore_existing_files=False)
#     symmp2.build_matrix(ignore_existing_files=False)
#     P1 = symmp1.get_matrix()
#     P2 = symmp2.get_matrix()
#     diff = P1*P2
#     diffs = sum(abs(c) for cc in diff.columns() for c in cc)
#     print("PPTest", n_vertices, n_loops,
#           n_hairs, n_ws, rep_ind1, rep_ind2, ", result: ", diffs)


def DDTest(n_vertices, n_loops, n_marked, n_hairs, even_edges, print_matrices=False, plot_bases=False):
    # tests whether differentials anti commute
    D1o = ForestedGraphComplex.ContractEdgesGO.generate_operator(
        n_vertices, n_loops, n_marked, n_hairs, even_edges)
    D2o = ForestedGraphComplex.ContractEdgesGO.generate_operator(
        n_vertices, n_loops, n_marked-1, n_hairs, even_edges)
    DD1o = ForestedGraphComplex.UnmarkEdgesGO.generate_operator(
        n_vertices, n_loops, n_marked, n_hairs, even_edges)
    DD2o = ForestedGraphComplex.UnmarkEdgesGO.generate_operator(
        n_vertices-1, n_loops, n_marked-1, n_hairs, even_edges)

    D1 = D1o.get_matrix()
    D2 = D2o.get_matrix()
    DD1 = DD1o.get_matrix()
    DD2 = DD2o.get_matrix()

    diff = DD2*D1+D2*DD1
    diffs = sum(abs(c) for cc in diff.columns() for c in cc)
    print("DDTest", n_vertices, n_loops,
          n_hairs, n_marked, even_edges, ", result: ", diffs)
    if diffs > 0 or print_matrices:
        print(D1)
        print(D2)
        print(DD1)
        print(DD2)
        print(diff)
    if plot_bases:
        D1o.domain.display_basis_plots()
        D1o.target.display_basis_plots()
        D2o.target.display_basis_plots()
        DD2o.target.display_basis_plots()


# def SumOneTest(n_vertices, n_loops, n_hairs, n_ws):
#     nparts = len(list(Partitions(n_hairs)))
#     Plist = []
#     for j in range(nparts):
#         symmp1 = WRHairyGraphComplex.SymmProjector.generate_operator(
#             n_vertices, n_loops, n_hairs, n_ws, j)
#         symmp1.build_matrix(ignore_existing_files=False)
#         P1 = symmp1.get_matrix()
#         P1 = P1.change_ring(
#             QQ) * SymmRepDimension(Partitions(n_hairs)[j]) / factorial(n_hairs)
#         Plist.append(P1)
#         # print(P1)

#     Psum = sum(Plist)
#     print(Psum)


# def getCohomDimP(n_vertices, n_loops, n_hairs, n_ws, rep_ind):
#     tt = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
#         n_vertices, n_loops, n_hairs, n_ws)
#     tu = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
#         n_vertices+1, n_loops, n_hairs, n_ws)
#     symmp1 = WRHairyGraphComplex.SymmProjector.generate_operator(
#         n_vertices, n_loops, n_hairs, n_ws, rep_ind)

#     D1 = tt.get_matrix()
#     D2 = tu.get_matrix()
#     # C = D2*D1
#     symmp1.build_matrix(ignore_existing_files=True)
#     P1 = symmp1.get_matrix()
#     print("matrices loaded")

#     D1P = D1*P1
#     D2P = P1*D2
#     print("computing ranks....")

#     diff = D1*D2
#     diffs = sum(abs(c) for cc in diff.columns() for c in cc)
#     print(diffs)

#     isocomp_dim = P1.rank()
#     r1 = D1P.rank()
#     r2 = D2P.rank()
#     print(isocomp_dim, r1, r2)
#     cohomdim = isocomp_dim - r1-r2
#     part = Partitions(n_hairs)[rep_ind]
#     rep_dim = symmetrica.charvalue(part, [1 for j in range(n_hairs)])
#     if cohomdim > 0:
#         print("Cohomology found:  w=", n_ws, ", h=", n_hairs, ", l=", n_loops, ", vertices=", n_vertices,
#               " (degree ", n_vertices+n_loops +
#               1, "), partition=", part,  ", invpartition=", part.conjugate(),
#               ", multiplicity=", cohomdim/rep_dim, ", cohomdim=", cohomdim)
#     return isocomp_dim - r1-r2


# def getCohomDimPAll(gvs):
#     for w in gvs.w_range:
#         for h in gvs.h_range:
#             for l in gvs.l_range:
#                 for v in gvs.v_range:
#                     D1 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
#                         v, l, h, w)
#                     D2 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(
#                         v+1, l, h, w)
#                     try:
#                         d = WRHairyGraphComplex.WRHairyGraphVS(
#                             v, l, h, w).get_dimension()
#                         r1 = D1.get_matrix_rank()
#                         r2 = D2.get_matrix_rank()
#                         if d-r1-r2 > 0:
#                             for rep_ind in range(len(Partitions(h))):
#                                 getCohomDimP(v, l, h, w, rep_ind)
#                     except:
#                         pass

maxl = 2
maxh = 1
for l in range(maxl+1):
    PFGC = ForestedGraphComplex.PreForestedGraphSumVS(
        range(2*maxl - 2 + maxh), range(l, l+1), range(2*maxl - 2+maxh-1), range(maxl-l+1+maxh))
    PFGC.build_basis(ignore_existing_files=True)
