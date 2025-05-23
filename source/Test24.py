# tests for Kneissler relations 
from sage.all import *
import KneisslerGC
from itertools import permutations

nloops = 9
k = nloops-1
even_edges = False

V0 = KneisslerGC.KneisslerGVS(nloops, 0, even_edges)
V1 = KneisslerGC.KneisslerGVS(nloops, 1, even_edges)
V2 = KneisslerGC.KneisslerGVS(nloops, 2, even_edges)
V3 = KneisslerGC.KneisslerGVS(nloops, 3, even_edges)

op0 = KneisslerGC.ContractEdgesGO.generate_operator(nloops, 0, even_edges)
op2 = KneisslerGC.ContractEdgesGO.generate_operator(nloops, 2, even_edges)
op3 = KneisslerGC.ContractEdgesGO.generate_operator(nloops, 3, even_edges)
D0 = op0.get_matrix()
D2 = op2.get_matrix()
D3 = op3.get_matrix()

bg60 = V0.get_basis_g6()
bdict0 = V0.get_g6_coordinates_dict() # barrel graphs
bdict1 = V1.get_g6_coordinates_dict() # 4-valent = tbarrel(=a) + xtbarrel(=b)
bdict2 = V2.get_g6_coordinates_dict() # all 3-valent = barrel(=a) + triangle(=b) + h(=c)
bdict3 = V3.get_g6_coordinates_dict() # complement of 0 in 2

# pieces of matrix D2:
barrel_idxs = [bdict2[s] for s in bg60]
tbarrelg6 = { V1.graph_to_canon_g6(G)[0] for G in KneisslerGC.all_tbarrel_graphs(k) }
xtbarrelg6 = { V1.graph_to_canon_g6(G)[0] for G in KneisslerGC.all_xtbarrel_graphs(k) }
tbarrel_idxs = [bdict1[s] for s in tbarrelg6 if s in bdict1 and not s in xtbarrelg6]
xtbarrel_idxs = [bdict1[s] for s in xtbarrelg6 if s in bdict1]
triangleg6 = { V3.graph_to_canon_g6(G)[0] for G in KneisslerGC.all_triangle_graphs(k) }
triangle_idxs = [bdict2[s] for s in triangleg6 if s in bdict2]# and not s in bg60]
hg6 = { V3.graph_to_canon_g6(G)[0] for G in KneisslerGC.all_hgraph_graphs(k) }
h_idxs = [bdict2[s] for s in hg6 if s in bdict2 and not s in bg60]

D2aa = D2[tbarrel_idxs, barrel_idxs]
D2ab = D2[tbarrel_idxs, triangle_idxs]
D2ac = D2[tbarrel_idxs, h_idxs]
D2ba = D2[xtbarrel_idxs, barrel_idxs]
D2bb = D2[xtbarrel_idxs, triangle_idxs]
D2bc = D2[xtbarrel_idxs, h_idxs]

print("Numbers of graphs in the different classes:")
print("barrel", len(barrel_idxs))
print("tbarrel", len(tbarrel_idxs))
print("xtbarrel", len(xtbarrel_idxs))
print("triangle", len(triangle_idxs))
print("h", len(h_idxs))

print("D2ab size", D2ab.nrows(), D2ab.ncols())
# for i in range(D2bb.ncols()):
#     print("D2bb #xts: ", len(D2bb.nonzero_positions_in_column(i)))

# for i in h_idxs:
    # print("h graph #nnz: ", len(D2.nonzero_positions_in_column(i)))
# for i in range(D2bc.ncols()):
#     print("D2bc #xts: ", len(D2bc.nonzero_positions_in_column(i)))

# for i in range(D2ab.ncols()):
#     print("D2ab #ts: ", len(D2ab.nonzero_positions_in_column(i)))

def leading_terms():
    ret = set()
    for p in permutations(range(k-1)):
        knrel = KneisslerGC.kneissler_rel(k, p, even_edges)
        knrel_g6 = [V0.graph_to_canon_g6(G)[0] for G,_ in knrel]
        knidxs = [bdict0[s] for s in knrel_g6 if s in bdict0]
        if len(knidxs) == 0:
            print("empty relation", knrel_g6)
            continue 
        # find lowest idx
        lowest_idx = min(knidxs)
        # check if lowest_idx appears only once in knidxs
        if knidxs.count(lowest_idx) == 1:
            ret.add(lowest_idx)
        else:
            print("leading term not unique: ", knrel_g6)
    return ret


ldterms = leading_terms()
print("leading terms len: ", len(ldterms), " of ", len(bg60))

# intersection checks
if (False):
# if (True):
    # check: intersection of triangle or h graphs with barrel?
    for s in triangleg6:
        if s in bg60:
            print("triangle in barrel", s)
    for s in hg6:
        if s in bg60:
            print("h in barrel", s)
    # check: intersection of triangle with h graphs
    for s in triangleg6:
        if s in hg6:
            print("triangle in h", s)
    # check: intersection of tbarrel with xtbarrel graphs
    for s in tbarrelg6:
        if s in xtbarrelg6:
            print("tbarrel in xtbarrel", s)

