import OrdinaryGraphComplex
import CHairyGraphComplex
from sage.all import *


def DSquareTestSingle(n_vertices, n_loops, even_edges, j_to_pick=-1, plot_basis=False):
    tt = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(
        n_vertices, n_loops, even_edges)
    tu = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(
        n_vertices-1, n_loops, even_edges)
    D1 = tt.get_matrix()
    D2 = tu.get_matrix()
    C = D2*D1

    print("D1")
    print(D1)
    print("D2")
    print(D2)
    print("D2*D1")
    print(C)

    ba0 = tt.domain.get_basis_g6()
    ba1 = tu.domain.get_basis_g6()
    ba2 = tu.target.get_basis_g6()

    if plot_basis:
        tt.domain.plot_all_graphs_to_file(skip_existing=False)
        tt.domain.display_basis_plots()
        tu.domain.plot_all_graphs_to_file(skip_existing=False)
        tu.domain.display_basis_plots()
        tu.target.plot_all_graphs_to_file(skip_existing=False)
        tu.target.display_basis_plots()

    if (j_to_pick < 0):
        for i in range(C.nrows()):
            for j in range(C.ncols()):
                if C[i, j] != 0:
                    print("Nonzero entry in product: ", i, j, C[i, j])
                    j_to_pick = j
        if j_to_pick < 0:
            print("success, squares to zero")
            return
        else:
            print("Does not square to zero, checking index ",
                  j_to_pick, " g6code ", ba0[j_to_pick])
    else:
        print("Checking index ",
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
                print(g6, " exists at index ", ba1.index(
                    g6), " v=", x, "sgn=", sgn)

    # compute D^2
    ww = [(HH, x*xx) for H, x in w for HH, xx in tu.operate_on(H)]
    wwd = {}
    for H, x in ww:
        g6, sgn = tu.target.graph_to_canon_g6(H)
        if g6 in wwd:
            wwd[g6] += x * sgn
        else:
            wwd[g6] = x * sgn
    print(wwd)
    nonzeroflag = false
    for g6, x in wwd.items():
        if x != 0:
            print("Nonzero entry: ", g6, x)
            nonzeroflag = true
    if not nonzeroflag:
        print("all entries zero, i.e., success.")


# DSquareTestSingle(8, 6, False, plot_basis=True)

OGC = OrdinaryGraphComplex.OrdinaryGC(
    range(17), range(9), True, {'contract'})

OGC.build_basis(ignore_existing_files=True)
OGC.build_matrix(ignore_existing_files=True)
OGC.compute_rank(sage='integer', ignore_existing_files=True)
# OGC.square_zero_test()
# OGC.test_pairwise_anti_commutativity()
OGC.export_cohomology_dim_for_web()

# DDD = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(
#     8, 6, False)

# w1 = DDD.operate_on(Graph("G@hU^_"))

# DDD2 = CHairyGraphComplex.ContractEdgesGO.generate_operator(
#     8, 6, 0, False)

# w2 = DDD2.operate_on(Graph("G@hU^_"))

# print([[DDD.target.graph_to_canon_g6(H), x] for H, x in w1])
# print([[DDD2.target.graph_to_canon_g6(H), x] for H, x in w2])
# print([[DDD.target.graph_to_canon_g6(H), x] for H, x in w2])
# print([[DDD2.target.graph_to_canon_g6(H), x] for H, x in w1])
# print([[H.graph6_string(), x] for H, x in w1])
# print([[H.graph6_string(), x] for H, x in w2])
# print(w1)
# print(w2)


# ttt = "FtYUW"
# # gtt = Graph(ttt)
# gtt, dummy = w2[1]
# print(gtt)
# # gtt = gtt.copy()
# VS1 = OrdinaryGraphComplex.OrdinaryGVS(7, 6, False)
# VS2 = CHairyGraphComplex.CHairyGraphVS(7, 6, 0, False)
# VS3 = DDD.target

# print(VS1.graph_to_canon_g6(gtt))
# print(VS2.graph_to_canon_g6(gtt))
# print(VS3.graph_to_canon_g6(gtt))
