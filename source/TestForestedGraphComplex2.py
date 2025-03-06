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
        for i in range(0, C.nrows()):
            for j in range(0, C.ncols()):
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


def loopg(nv):
    g = Graph(4*nv)
    for i in range(nv):
        # marked edges
        g.add_edge(2*i, 2*i+1)
        # unmarked edges parallel to marked edge
        g.add_edge(2*nv+i, 2*i)
        g.add_edge(2*nv+i, 2*i+1)
        # other unmarked edge
        g.add_edge(3*nv+i, 2*i+1)
        g.add_edge(3*nv+i, (2*i+2) % (2*nv))
    return g

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


looplen = 5

D1o = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
    looplen+1, looplen, 0, True)
D2o = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
    looplen+1, looplen+1, 0, True)
vs = ForestedGraphComplex.ForestedGVS(looplen * 2, looplen+1, looplen, 0, True)

# D1o.domain.plot_all_graphs_to_file(skip_existing=False)
# D1o.domain.display_basis_plots()

D1 = D1o.get_matrix()
D2 = D2o.get_matrix()

ba1 = D1o.domain.get_basis_g6()

print("D1: ", D1.nrows(), D1.ncols())
print("D2: ", D2.nrows(), D2.ncols())
# print(ba1[37])

gg = loopg(looplen)
ggs = gg.canonical_label(partition=vs.get_partition()).graph6_string()
# print(ggs)
ggi = ba1.index(ggs)

print("Index of loop graph in basis: ", ggi)
# gg.show()
# gg = loopg(4)
# gg.show()

ppp = random_prime(2222222222, False, 222222222)
print("Prime: ", ppp)

vloop = matrix(len(ba1), 1)
vloop[ggi] = 1

D2f = D2.change_ring(RDF)
vloopf = vloop.change_ring(RDF)
print("Done1")

xxx = D2f.solve_right(vloopf)
print("Done2")
print((vloopf-D2f*xxx).norm())

# print(D1 * vloop)
# # print("v:", vloop)

# y = matrix(D2.ncols()+1, 1)
# y[D2.ncols()] = 1
# D2tx = D2.transpose().stack(vloop.transpose())
# print("D2tx: ", D2tx.nrows(), D2tx.ncols())
# z = D2tx.solve_right(y)
# print(z)

#sol = D2.solve_right(vloop)

# print(D2*sol)


# print(D1)
# print(D2)

# tt = ForestedGraphComplex.PreForestedGVS(4, 4, 1, 0)
# # tt = ForestedGraphComplex.PreForestedGVS(4, 3, 3, 2)
# tt.build_basis(ignore_existing_files=True)
# tt.plot_all_graphs_to_file(skip_existing=False)
# tt.display_basis_plots()

# ttt = ForestedGraphComplex.ForestedGVS(4, 4, 1, 0, False)
# tt = ForestedGraphComplex.PreForestedGVS(4, 3, 3, 2)
# ttt.build_basis(ignore_existing_files=True)
# ttt.plot_all_graphs_to_file(skip_existing=False)
# ttt.display_basis_plots()

# GG = Graph("I@XbCaOc?")
# autom_list = GG.automorphism_group(partition=ttt.get_partition()).gens()
# print(ttt._has_odd_automorphisms(GG, autom_list))

# print(diff)
# print(P)

# v_range = range(1, 8)
# l_range = range(0, 7)
# h_range = range(1, 8)
# edges_types = [True, False]
# hairs_types = [True, False]

# tt = WRHairyGraphComplex.WRHairyGraphVS(2, 2, 2, 2)
# tt.build_basis()
# # tt.plot_all_graphs_to_file()
# tt.display_basis_plots()
# tt2 = WRHairyGraphComplex.WRHairyGraphVS(1, 2, 2, 2)
# tt2.build_basis()
# # tt.plot_all_graphs_to_file()
# tt2.display_basis_plots()

# tt = WHairyGraphComplex.WHairyGraphVS(4,3,0,1)
# print(tt.is_valid())
# tt.build_basis(ignore_existing_files=True)

# tt.plot_all_graphs_to_file()
# tt.display_basis_plots()


# WGC = WHairyGraphComplex.WHairyGC(range(0,11), range(5,7), range(0,4), range(2,3) , ['contract'])
# WGC = WRHairyGraphComplex.WRHairyGC(range(0, 14), range( 0, 6), range(2, 3), range(2, 3), ['contract'])
# WGC = WRHairyGraphComplex.WRHairyGC(range(0,10), range(0,2), range(4,7), range(1,2) , ['contract'])
# WGC = WHairyGraphComplex.WHairyGC(range(0,8), range(0,6), range(1,3), range(2,3) , ['contract'])

# maxl = 6
# for l in range(maxl+1):
#     PFGC = ForestedGraphComplex.PreForestedGraphSumVS(
#         range(0, 2*maxl - 1), range(l, l+1), range(0, 2*maxl - 2), range(0, maxl-l+1))
#     PFGC.build_basis(ignore_existing_files=False)

# PFGC = ForestedGraphComplex.PreForestedGraphSumVS(
#     range(0, 11), range(0, 7), range(0, 10), range(0, 1))
# PFGC.build_basis(ignore_existing_files=False)
# PFGC = ForestedGraphComplex.PreForestedGraphSumVS(
#     range(0, 8), range(0, 4), range(0, 7), range(1, 2))
# PFGC.build_basis(ignore_existing_files=True)
# PFGC = ForestedGraphComplex.PreForestedGraphSumVS(
#     range(0, 8), range(0, 3), range(0, 7), range(2, 4))
# PFGC.build_basis(ignore_existing_files=True)
# PFGC = ForestedGraphComplex.PreForestedGraphSumVS(
#     range(0, 8), range(0, 1), range(0, 7), range(4, 5))
# PFGC.build_basis(ignore_existing_files=True)

# FGC = ForestedGraphComplex.ForestedGC(
#     range(0, 11), range(0, 7), range(0, 10), range(0, 1), True, {'contract', 'unmark'})
# FGC.build_basis(ignore_existing_files=False)
# FGC.build_matrix(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=False)
# FGC.square_zero_test()

# DDTest(4, 3, 2, 0, False, plot_bases=False)

# for v in range(2, 8):
#     for g in range(5):
#         for m in range(2, 7):
#             DDTest(v, g, m, 0, False)

# FBGC = ForestedGraphComplex.ForestedContractUnmarkBiGC(
#     range(0, 7), range(0, 10), range(0, 1), True)
# FBGC.build_basis(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=False)
# FBGC.build_matrix(progress_bar=False, info_tracker=False,
#                   ignore_existing_files=False)
# # FBGC.square_zero_test()

# FBGC.compute_rank(ignore_existing_files=False, sage="integer")


# FBGC.print_dim_and_eulerchar()
# FBGC.print_cohomology_dim()

# l = 3
# m = 0
# h = 0
# even_edges = False
# D1 = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
#     l, m, h, even_edges)
# D2 = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
#     l, m+1, h, even_edges)
# d = ForestedGraphComplex.ForestedDegSlice(
#     l, m, h, even_edges).get_dimension()
# r1 = D1.get_matrix_rank()
# r2 = D2.get_matrix_rank()
# print(d-r1-r2)


# DSquareTestSingleUnmark(4, 4, 2, 0, False, plot_basis=True)


# WGC = ForestedGraphComplex.ForestedGC(range(0, 14), range(
#     1, 2), range(4, 6), range(1, 2), ['contract'])

# WGC = WRHairyGraphComplex.WRHairyGC(range(0, 14), range(
#     0, 1), range(4, 5), range(1, 2), ['contract'])

# WGC.build_basis(progress_bar=False, info_tracker=False,
#                 ignore_existing_files=True)
# WGC.build_matrix(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=True)

# WGC.build_basis(progress_bar=False, info_tracker=False,
#                 ignore_existing_files=True)
# WGC.build_matrix(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=True)

# # WGC.square_zero_test()

# # WGC.compute_rank(ignore_existing_files=True, sage="mod")
# WGC.compute_rank(ignore_existing_files=True, sage="integer")
# # WGC.plot_cohomology_dim(to_html=True)
# # Euler char


# v = 2
# l = 0
# h = 4
# w = 1

# SumOneTest(v, l, h, w)
# for j in range(len(Partitions(h))):
#     PSquareTest(v, l, h, w, j)
#     PDTest(v, l, h, w, j)


# for j in range(len(Partitions(h))):
#     print(getCohomDimP(v, l, h, w, j))

# PSquareTest(4, 3, 2, 2, 0)
# PSquareTest(3, 3, 2, 2, 0)

# PPTest(4, 3, 2, 2, 0, 1)
# PDTest(4, 3, 2, 2, 0)

# PDTest(2, 2, 2, 2, 0, print_matrices=True)
# PDTest(3, 2, 2, 2, 0)
# PDTest(4, 2, 2, 2, 0)
# PDTest(4, 3, 2, 2, 0)
# PDTest(3, 3, 2, 2, 0)

# print(getCohomDimP(6, 5, 2, 2, 0))
# print(getCohomDimP(8, 5, 2, 2, 1))

# print(getCohomDimP(6, 5, 2, 2, 1))

# SumOneTest(4, 2, 3, 1)
# PSquareTest(4, 2, 3, 1, 0)
# PSquareTest(4, 2, 3, 1, 1)
# PSquareTest(4, 2, 3, 1, 2)
# print(getCohomDimP(2, 1, 5, 2, 1))
# print(getCohomDimP(4, 1, 5, 2, 2))
# print(getCohomDimP(4, 1, 5, 2, 3))
# print(getCohomDimP(4, 1, 5, 2, 4))
# print(getCohomDimP(4, 1, 5, 2, 5))

# getCohomDimPAll(WGC)

# print(getCohomDimP(3, 0, 6, 2, 0))
# print(getCohomDimP(3, 0, 6, 2, 1))
# print(getCohomDimP(3, 0, 6, 2, 2))
# print(getCohomDimP(3, 0, 6, 2, 3))
# print(getCohomDimP(3, 0, 6, 2, 4))
# print(getCohomDimP(3, 0, 6, 2, 5))
# print(getCohomDimP(3, 0, 6, 2, 6))
# print(getCohomDimP(3, 0, 6, 2, 7))
# print(getCohomDimP(3, 0, 6, 2, 8))
# print(getCohomDimP(3, 0, 6, 2, 9))
# print(getCohomDimP(3, 0, 6, 2, 10))

# print(getCohomDimP(3, 1, 4, 2, 0))
# print(getCohomDimP(3, 1, 4, 2, 1))
# print(getCohomDimP(3, 1, 4, 2, 2))
# print(getCohomDimP(3, 1, 4, 2, 3))
# print(getCohomDimP(3, 1, 4, 2, 4))

# print(getCohomDimP(2, 0, 5, 2, 0))
# print(getCohomDimP(2, 0, 5, 2, 1))
# print(getCohomDimP(2, 0, 5, 2, 2))
# print(getCohomDimP(2, 0, 5, 2, 3))
# print(getCohomDimP(2, 0, 5, 2, 4))
# print(getCohomDimP(2, 0, 5, 2, 5))
# print(getCohomDimP(2, 0, 5, 2, 6))

# print(getCohomDimP(4, 1, 5, 2, 0))
# print(getCohomDimP(4, 1, 5, 2, 1))
# print(getCohomDimP(4, 1, 5, 2, 2))
# print(getCohomDimP(4, 1, 5, 2, 3))
# print(getCohomDimP(4, 1, 5, 2, 4))
# print(getCohomDimP(4, 1, 5, 2, 5))
# print(getCohomDimP(4, 1, 5, 2, 6))


# print(getCohomDimP(7, 4, 3, 1, 0))
# print(getCohomDimP(7, 4, 3, 1, 1))
# print(getCohomDimP(7, 4, 3, 1, 2))
# print(getCohomDimP(7, 4, 3, 2, 0))
# print(getCohomDimP(7, 4, 3, 2, 1))
# print(getCohomDimP(7, 4, 3, 2, 2))
# print(getCohomDimP(8, 5, 2, 2, 1))
# print(getCohomDimP(4, 3, 2, 1, 2))
# print(getCohomDimP(3, 1, 3, 1, 2))

# symmp = WRHairyGraphComplex.SymmProjector.generate_operator(4, 3, 2, 2, 0)
# symmp.build_matrix(ignore_existing_files=True)
# P = symmp.get_matrix()
# # diff = P*P-2*P  # should be zero
# # print(diff)
# # print(P)
# symmp2 = WRHairyGraphComplex.SymmProjector.generate_operator(3, 3, 2, 2, 0)
# symmp2.build_matrix(ignore_existing_files=True)
# P2 = symmp2.get_matrix()

# tt = WRHairyGraphComplex.ContractEdgesGO.generate_operator(4, 3, 2, 2)
# D1 = tt.get_matrix()
# diff = P2*D1-D1*P  # should be zero
# print(diff)

# WGC.plot_info()

# tt = WRHairyGraphComplex.WRHairyGraphVS(3,3,0,2)
# tt.build_basis(ignore_existing_files=True)
# tt.plot_all_graphs_to_file(skip_existing=False)
# tt.display_basis_plots()

# DSquareTestSingle(4,3,0,2, plot_basis=False)


# HG = WHairyGraphComplex.WHairyGraphVS(7,4,2,2)
# huntfor = "JSSqW@_?c??"
# huntforG = Graph(huntfor)

# ba = HG.get_basis_g6()
# print("In basis?: ", huntfor in ba)
# autom_list = huntforG.automorphism_group(partition=HG.get_partition()).gens()
# print("Odd automs?: ", HG._has_odd_automorphisms(huntforG,autom_list))
# print("auts", autom_list)
# print(huntforG)
# for p in autom_list:
#     print(list(p.tuple()))

# genlist=HG.get_generating_graphs()
# genlistg6 =[HG.graph_to_canon_g6(G)[0] for G in genlist]
# print("In gen list?: ", huntfor in genlistg6)
# testg6=HG.graph_to_canon_g6(huntforG)[0]
# print(testg6, testg6 in genlistg6)

# all_perm = [ list(range(0,HG.n_vertices+2)) + list(p) for p in itertools.permutations(range(HG.n_vertices+2, HG.n_vertices+HG.n_hairs+2)) ]
# print(all_perm)

# huntforG_h = huntforG.relabel({7:8, 8:7}, inplace=false)
# huntforG_h.add_edge(7,8)
# hairy_part= [list(range(HG.n_vertices+1)), list(range(HG.n_vertices+1, HG.n_vertices+HG.n_hairs+2))]
# huntforG_hg6 = huntforG_h.canonical_label(partition=hairy_part)

# somelist = HG._hairy_to_w_hair_2w(huntforG_h)
# somelistg6 =[HG.graph_to_canon_g6(G)[0] for G in somelist]
# print(huntfor in somelistg6)
# print(somelistg6)

# HG.build_basis(ignore_existing_files=True)

# check_graphs_vs_basis(tu.target, wwd)

# tt.domain.plot_all_graphs_to_file(skip_existing=False)
# tt.domain.display_basis_plots()
# tu.domain.plot_all_graphs_to_file(skip_existing=False)
# tu.domain.display_basis_plots()
# tu.target.plot_all_graphs_to_file(skip_existing=False)
# tu.target.display_basis_plots()


# g = Graph("EZq?")
# res = tt.operate_on(g)
# print(res)
# for G,v in res:
#     g6, s = tt.target.graph_to_canon_g6(G)
#     print(g6, v*s)


# class BasisTest(TestGraphComplex.BasisTest):
#     def setUp(self):
#         self.vs_list = [HairyGraphComplex.HairyGraphVS(v, l, h, even_edges, even_hairs) for (v, l, h, even_edges, even_hairs)
#                         in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


# class OperatorTest(TestGraphComplex.OperatorTest):
#     def setUp(self):
#         self.op_list = [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, even_edges, even_hairs) for
#                         (v, l, h, even_edges, even_hairs) in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


# class GraphComplexTest(TestGraphComplex.GraphComplexTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]


# class CohomologyTest(TestGraphComplex.CohomologyTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]


# class SquareZeroTest(TestGraphComplex.SquareZeroTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]

# class AntiCommutativityTest(TestGraphComplex.AntiCommutativityTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, False, False, ['contract', 'et1h'])]

# def suite():
#     suite = unittest.TestSuite()
#     suite.addTest(BasisTest('test_basis_functionality'))
#     suite.addTest(BasisTest('test_compare_ref_basis'))
#     suite.addTest(OperatorTest('test_operator_functionality'))
#     suite.addTest(OperatorTest('test_compare_ref_op_matrix'))
#     suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
#     suite.addTest(CohomologyTest('test_cohomology_functionality'))
#     suite.addTest(SquareZeroTest('test_square_zero'))
#     suite.addTest(AntiCommutativityTest('test_anti_commutativity'))
#     return suite


# if __name__ == '__main__':
#     print("\n#######################################\n" + "----- Start test suite for hairy graph complex -----")
#     runner = unittest.TextTestRunner()
#     runner.run(suite())
