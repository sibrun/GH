import OrdinaryVariants
import OrdinaryGraphComplex
from sage.all import *


def get_filtered_matrix_rank(A, F1, F2, f):
    """Computes the rank of the part of matrix A consisting of those rows where list F1 takes value f and those columsn where F2 takes value f

    :param A: _description_
    :type A: _type_
    :param F1: _description_
    :type F1: _type_
    :param F2: _description_
    :type F2: _type_
    :param f: _description_
    :type f: _type_
    """
    # find list of indices where F1 is f
    idx1 = [i for i, x in enumerate(F1) if x == f]
    # find list of indices where F2 is f
    idx2 = [i for i, x in enumerate(F2) if x == f]
    if len(idx1) == 0 or len(idx2) == 0:
        return 0
    # extract the submatrix
    B = A.matrix_from_rows_and_columns(idx1, idx2)
    # compute the rank
    return B.rank()

def get_filtered_op_rank(D, F1, F2, f):
    if not D.is_valid():
        return 0
    
    if D.get_domain().get_dimension() == 0 or D.get_target().get_dimension() == 0:
        return 0
    
    return get_filtered_matrix_rank(D.get_matrix(), F1, F2, f)

def get_vs_filter_list(V, filter_func):
    if not V.is_valid():
        return []
    else:
        return [filter_func(G) for G in V.get_basis()]

def print_filtered_cohomology(vs_class, op_class, l, even_edges, filter_func, flist, ignore_ex = False):
    # make list of vector spaces 
    # allvs = [OrdinaryVariants.OrdinaryGVSFull(v, l, even_edges) for v in range(2, 21)]
    allvs = [vs_class(v, l, even_edges) for v in range(2, 21)]
    for V in allvs:
        V.build_basis(ignore_existing_files=ignore_ex)
    # make list of operators
    allops = [op_class.generate_operator(v, l, even_edges) for v in range(3, 21)]
    # allops = [OrdinaryVariants.ContractEdgesGOFull.generate_operator(v, l, even_edges) for v in range(3, 21)]
    for D in allops:
        D.build_matrix(ignore_existing_files=ignore_ex)
    # create all filter lists
    filter_list = [get_vs_filter_list(V, filter_func) for V in allvs]
    # print(filter_list)
    # zipped filter pairs for operators
    filter_zip = list(zip(filter_list[:-1], filter_list[1:]))

    for f in flist:
        print("Filter ", f)
        # compute matrix ranks
        ranks = [get_filtered_op_rank(D, F1, F2, f) for D, (F1, F2) in zip(allops, filter_zip)]
        # compute cohomology dimensions
        fil_counts = [fil.count(f) for fil in filter_list]
        print(fil_counts)
        rank_zip = list(zip(fil_counts[1:-1],ranks[:-1], ranks[1:]))
        # print(rank_zip)
        cohom = { i+3: fcount - r1 - r2 for i, (fcount, r1, r2) in enumerate(rank_zip)}
        print(cohom)

def vertex_connectivity_filter(G):
    return G.vertex_connectivity()

def triangle_filter(G):
    # number of triangles in the graph
    return len(list(sage.graphs.cliquer.all_cliques(G, 3,3)))

def atoms_filter(G):
    atoms, cliques = G.atoms_and_clique_separators()
    return len(atoms)

def highvalence_filter(G):
    return sum( (1 if len(G[v])>4 else 0) for v in G.vertices(sort=True) ) 

def maxvalence_filter(G):
    return max( len(G[v])  for v in G.vertices(sort=True) ) 


# print_filtered_cohomology(6, False, triangle_filter, [0, 1,2,3,4, 5], ignore_ex=False)
# print_filtered_cohomology(OrdinaryGraphComplex.OrdinaryGVS, OrdinaryGraphComplex.ContractEdgesGO, 8, False, atoms_filter, [0, 1,2,3,4, 5], ignore_ex=True)
# print_filtered_cohomology(OrdinaryVariants.OrdinaryGVSTriconnected, OrdinaryVariants.ContractEdgesGOTriconnected, 8, False, atoms_filter, [0, 1,2,3,4, 5], ignore_ex=True)
# print_filtered_cohomology(OrdinaryVariants.OrdinaryGVSTriconnected, OrdinaryVariants.ContractEdgesGOTriconnected, 8, False, triangle_filter, [0, 1,2,3,4, 5], ignore_ex=False)
# print_filtered_cohomology(OrdinaryVariants.OrdinaryGVSTriconnected, OrdinaryVariants.ContractEdgesGOTriconnected, 8, False, highvalence_filter, [0, 1,2,3,4, 5], ignore_ex=False)
print_filtered_cohomology(OrdinaryVariants.OrdinaryGVSFull, OrdinaryVariants.ContractEdgesGOFull, 8, False, maxvalence_filter, [3,4, 5,6,7,8], ignore_ex=False)
# print_filtered_cohomology(8, False, vertex_connectivity_filter, [1,2,3,4])



# test triangle counts 
# print(get_vs_filter_list(OrdinaryGraphComplex.OrdinaryGVS(10, 6, False), triangle_filter))