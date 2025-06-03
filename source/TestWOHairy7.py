# Testing new version of WOHary
from sage.all import *
import WOHairyGC
import WOHairyGraphComplex2



def fuse_epsilon(G, n_ws, n_epsilons, n_numbered):
    """Translate a graph from Skipness representation to mine
    Skipness : [internal] [numbered] [omegas] [epsilons]
    My       : [internal] [single epsilon] [omegas] [numbered]
    """
    n_internal = G.order() - n_ws - n_epsilons - n_numbered
    n_edges = G.size()
    G1 = G.copy()
    if n_epsilon == 0:
        # add one epsilon vertex
        G1.add_vertex()
    else:
        G1.merge_vertices(range(n_internal+n_numbered+n_ws, n_internal+n_numbered+n_ws+n_epsilons))
    
    G1.relabel(list(range(n_internal)) +[i+n_internal+n_ws+1 for i in range(n_numbered)] + [i+n_internal+1 for i in range(n_ws)]  + [n_internal])
    return G1


def compare_bases(n_internal, n_loops, n_numbered, n_ws):
    n_edges = n_loops + n_internal
    degree = n_edges + 22 - n_ws
    VP = WOHairyGC.WOHairyGVS(n_loops+1, n_numbered, n_ws, degree)
    
    VT = WOHairyGraphComplex2.WOHairyGraphVS(n_internal, n_loops, n_numbered, n_ws)
    lst = VT.get_prerequisites()
    for U in lst:
        for UU in U.get_chairy_prerequisites():
            UU.build_basis()
        for UU in U.get_self_prerequisites():
            UU.build_basis()
        U.build_basis()
    
    VP.build_basis()
    VT.build_basis()
    print("VP: ", VP.get_dimension())
    print("VT: ", VT.get_dimension())

# compare_bases(3, 6, 5, 11)

VP = WOHairyGC.WOHairyGVS(9, 1, 11, 26)
VP.build_basis()
# print(VP.get_basis_file_path())
# for s in VP.get_basis_g6():
#     print(s)

# VP.display_basis_plots()

VPC = WOHairyGC.WOHairyComponentGVS(1, 0, 1, 2, 1)
VPC.build_basis()
VPC.display_basis_plots()