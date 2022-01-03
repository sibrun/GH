import OrdinaryGraphComplex
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
        tt.domain.display_basis_plots()
        tu.domain.display_basis_plots()
        tu.target.display_basis_plots()

    if (j_to_pick < 0):
        for i in range(0, C.nrows()):
            for j in range(0, C.ncols()):
                if C[i, j] != 0:
                    print("Nonzero entry in product: ", i, j, C[i, j])
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


DSquareTestSingle(8, 6, True, plot_basis=True)
