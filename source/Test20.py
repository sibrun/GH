# Test of Kneisslers method 
from sage.all import *
from itertools import permutations
import OrdinaryGraphComplex


def barrel_graph(k, p):
    # generates the barrel graph of 2k vertices, with p a permutation of the numbers 0,..,k-2
    G = Graph(2*k)
    # generate rims of barrel
    for j in range(k):
        G.add_edge(j, (j+1) % k)
        G.add_edge(k+j, k+(j+1)%k)
    # generate spokes
    G.add_edge(k-1, 2*k-1)
    for i,j in enumerate(p):
        G.add_edge(i, k+j)

    return G

def all_barrel_graphs(k):
    # generates all barrel graphs of 2k vertices
    for p in permutations(range(k-1)):
        yield barrel_graph(k, p)

def barrel_indices(k, even_edges):
    # finds the indices of the barrel graphs in the basis of the appropriate vector space
    V = OrdinaryGraphComplex.OrdinaryGVS(2*k, k+1, even_edges)
    B = V.get_g6_coordinates_dict()
    ret = set()
    for G in all_barrel_graphs(k):
        g6,_ = V.graph_to_canon_g6(G)
        if g6 in B:
            # append B[g6] to ret
            ret.add(B[g6])
    return ret

def get_neighbors(D, ids):
    # D is a matrix and ids a set of column indices.
    # returns a set of row indices containing those rows that share a nonzero with a column in ids
    ret = set()
    for i in ids:
        for j in D.nonzero_positions_in_column(i):
            ret.add(j)
    return ret

def get_neighbors2(D, ids):
    # D is a matrix and ids a set of row indices.
    # returns a set of row indices containing those rows that have a non-zero in a column in which ids aso have a nonzero
    return get_neighbors(D, get_neighbors(D.transpose(), ids))

def tbarrel_graph(k,p):
    # barrel graph with one 4-valent vertex (total 2k-1 vertices)
    G = barrel_graph(k, p)
    G.merge_vertices([k-2,k-1])
    G.relabel(list(range(0, G.order())), inplace=True)
    return G

def all_tbarrel_graphs(k):
    # generates all barrel graphs of 2k vertices
    for p in permutations(range(k-1)):
        yield tbarrel_graph(k, p)

def tbarrel_indices(k, even_edges):
    # finds the indices of the barrel graphs in the basis of the appropriate vector space
    V = OrdinaryGraphComplex.OrdinaryGVS(2*k-1, k+1, even_edges)
    V.build_basis()
    B = V.get_g6_coordinates_dict()
    ret = set()
    for G in all_tbarrel_graphs(k):
        g6,_ = V.graph_to_canon_g6(G)
        if g6 in B:
            # append B[g6] to ret
            ret.add(B[g6])
    return ret

def check_kneissler(k, even_edges):
    G = OrdinaryGraphComplex.OrdinaryGVS(2*k, k+1, even_edges)
    G.build_basis()
    op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(2*k, k+1, even_edges)
    op.build_matrix()
    op.compute_rank(sage="integer")
    D = op.get_matrix()
    m = D.nrows()
    n = D.ncols() # number of trivalent graphs
    Dt = D.transpose()
    bids = barrel_indices(k, even_edges)
    tids = get_neighbors(D, bids)
    #tids = get_neighbors2(D, tids)
    #tids = tbarrel_indices(k, even_edges)
    tids = get_neighbors2(D, tids)
    Dt2 = Dt[:,list(tids)]
    not_bids = set(range(Dt2.nrows())) - bids
    B = Dt2[list(not_bids), :]

    r0 = op.get_matrix_rank()
    r1 = Dt2.rank()
    r2 = B.rank()
    print("k=", k, "ee=", even_edges, ": Actual cohomdim: ", n-r0, " estimated cohomdim: ", len(bids) - r1 + r2, "(barrels: ", len(bids), " tbarrels:", len(tids), " non-barrels: ", len(not_bids), " non-tbarrels:", m-len(tids),")")
    
kk = 4
G = barrel_graph(kk, [0,2,1])
VV = OrdinaryGraphComplex.OrdinaryGVS(2*kk, kk+1, True)
g6,_ = VV.graph_to_canon_g6(G)
print("G6: ", g6)
print("Basis of GVS: ", VV.get_basis_g6())


# check_kneissler(3, True)
# check_kneissler(4, True)
# check_kneissler(5, True)
# check_kneissler(6, True)
# check_kneissler(7, True)
# check_kneissler(8, True)

check_kneissler(3, False)
check_kneissler(4, False)
check_kneissler(5, False)
check_kneissler(6, False)
check_kneissler(7, False)
check_kneissler(8, False)
check_kneissler(9, False)