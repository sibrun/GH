import unittest
import itertools
import logging
import Log
import TestGraphComplex
from WOHairyGraphComplex import *
from sage.all import *


def DSquareTestSingle(D1_go, D2_go, j_to_pick=-1, plot_basis=False):
    tt = D1_go
    tu = D2_go
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
        # tt.domain.display_basis_plots()
        # tu.domain.display_basis_plots()
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
            if not g6 in ba1:
                print(g6, " not found in basis ", " v=", x, " sgn=", sgn)
            else:
                print(g6, " exists at index ", ba1.index(
                    g6), " v=", x, " sgn=", sgn)

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


# GC = WOHairyGC(range(7), range(5), range(3), range(3), ['contract'])
# GC = WOHairyGC(range(5), range(3), range(3), range(3), ['epstoomega'])
GC = WOHairyGC(range(6), range(3,4), range(
    4), range(5), ['epstoomega', 'contract'])


GC.build_basis(ignore_existing_files=False)
# GC.build_basis(ignore_existing_files=True)

GC.build_matrix(ignore_existing_files=False)
# GC.build_matrix(ignore_existing_files=True)

GC.compute_rank(sage="integer")

# VS1 = WOHairyGraphVS(1, 4, 1, 3)
# print(VS1.is_valid())
# VS1.build_basis(ignore_existing_files=True)
# VS1.plot_all_graphs_to_file(skip_existing=False)
# VS1.display_basis_plots()

# GC.square_zero_test()

# GC.test_pairwise_anti_commutativity()

GC.print_dim_and_eulerchar()
GC.print_cohomology_dim()



go1 = ContractEdgesGO.generate_operator(5,3,2,1)
go2 = ContractEdgesGO.generate_operator(4,3,2,1)
# go2.target.plot_all_graphs_to_file(skip_existing=False)
# DSquareTestSingle(go1, go2, plot_basis=True)

D = go1.get_matrix()
DD = go2.get_matrix()

print(DD.right_kernel())
print(D.transpose().image())

go2.domain.display_basis_plots()
# G1 = Graph("EkQ?")
# ret = go1.operate_on(G1)
# for g, x in ret:
#     g6, sgn = go1.target.graph_to_canon_g6(g)
#     print(g.graph6_string(), g6, x, sgn)

# G2 = Graph(4)
# G2.add_edge((0,1))
# G2.add_edge((0,2))
# G2.add_edge((0,3))
# go3 = ContractEdgesGO.generate_operator(1,1,1,1)
# ret = go3.operate_on(G2)
# for g, x in ret:
#     print(g.graph6_string(), x)
