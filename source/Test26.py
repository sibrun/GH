# compute cohomology of forested 2 hair with spine graph complex

from sage.all import *
import ForestedGraphComplex

g = 2
nv = 2*g
m = 2

V1 = ForestedGraphComplex.ForestedGVS(nv, g, m, 2, False)
V1.build_basis()
b1 = V1.get_basis()
opc1 = ForestedGraphComplex.ContractEdgesGO.generate_operator(nv, g, m, 2, False)
opc1.build_matrix()
opuu = ForestedGraphComplex.UnmarkEdgesGO.generate_operator(nv, g, m, 2, False)
opuu.build_matrix()
Dc1 = opc1.get_matrix()
# print(Dc1)


opc2 = ForestedGraphComplex.ContractEdgesGO.generate_operator(nv, g, m+1, 2, False)
opc2.build_matrix()
V2 = opc2.domain
Dc2 = opc2.get_matrix()
b2 = V2.get_basis()
# print(Dc2)

def deleteUnmarkedEdges(G, V):
    G = G.copy()
    n = V.n_vertices
    m = V.n_unmarked_edges
    for i in range(n,n+m):
        for j in G.neighbors(i):
            G.delete_edge(i,j)
    return G

def are_hairs_connected(G,V):
    G = deleteUnmarkedEdges(G, V)
    n = V.n_vertices
    m = V.n_unmarked_edges
    h = V.n_hairs
    for i in range(n+m, n+m+h):
        for j in range(i+1, n+m+h):
            # print(i,j,G.distance(i,j))
            if G.distance(i,j) == Infinity:
                return False
    return True


def is_admissible(G,V):
    return are_hairs_connected(G, V)
    

def basis_filter(V):
    return [is_admissible(g, V) for g in V.get_basis()]

def filter_cols(M, fil):
    cols = [i for i, ok in enumerate(fil) if ok]
    return M[:, cols]
def filter_rows(M, fil):
    rows = [i for i, ok in enumerate(fil) if ok]
    return M[rows, :]
def filter_both(M, fil_rows, fil_cols):
    rows = [i for i, ok in enumerate(fil_rows) if ok]
    cols = [i for i, ok in enumerate(fil_cols) if ok]
    return M[rows, cols]

def filter_from_sumvs(V, filterfunc):
    return [filterfunc(g, VV) for VV in V.get_vs_list() for g in VV.get_basis() ]

def cohom_formatted_forested_top(D1, D2, Dc2, filterfunc):
    # D1.get_domain().build_basis()
    # D2.get_domain().build_basis()
    # Dc2.get_domain().build_basis()
    # D1.get_target().build_basis()
    # D2.get_target().build_basis()
    # Dc2.get_target().build_basis()
    # D1.build_matrix()
    # D2.build_matrix()
    # Dc2.build_matrix()

    vs = D1.get_domain()
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    
    vs_fil = [filterfunc(g, vs) for g in vs.get_basis()]
    out1_fil = filter_from_sumvs(D1.get_target(), filterfunc)

    # d = vs.get_dimension()
    d = sum( 1 for ok in vs_fil if ok)
    # print("dim", d)
    # print(vs_fil)

    r1 = 0
    r2 = 0
    rc2 = 0

    vs2 = D2.get_domain()
    vs2_fil = [filterfunc(g, vs2) for g in vs2.get_basis()]
    out2_fil = filter_from_sumvs(D2.get_target(), filterfunc)
    d2 = sum( 1 for ok in vs2_fil if ok)
    
    if Dc2.is_valid():
  
        Ac2 = Dc2.get_matrix()
        Ac2_fil = filter_cols(Ac2, vs2_fil)

        rc2 = Ac2_fil.rank()
        # print("Ac2", Ac2.nrows(), Ac2.ncols(), rc2)
        # print("Ac2_fil", Ac2_fil.nrows(), Ac2_fil.ncols(), rc2)
        

    if D1.is_valid():
        A1 = D1.get_matrix()
        A1_fil = filter_both(A1, out1_fil, vs_fil)
        # A1_fil = filter_cols(A1, vs_fil)
        r1 = A1_fil.rank()
        # print("A1", A1.nrows(), A1.ncols(), r1)
        # print("A1_fil", A1_fil.nrows(), A1_fil.ncols(), r1)

    if D2.is_valid():
        A2 = D2.get_matrix()
        A2_fil = filter_both(A2, out2_fil, vs2_fil)
        r2 = A2_fil.rank()
        # print("A2", A2.nrows(), A2.ncols(), r2)
        # print("A2_fil", A2_fil.nrows(), A2_fil.ncols(), r2)

    # exact or not?
    r_str = "" #f"({d}+{rc2}-{r1}-{r2},{d},{d2})"

    # iso string
    cohomdim = d+rc2-r1-r2

    return str(cohomdim) + r_str 


def create_forested_top_cohom_table(l_range, m_range, even_edges, filterfunc):
    s = ""
    h=2
    for l in l_range:
        print("l =", l)
        print(  [cohom_formatted_forested_top(
                    ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
                        l, m, h, even_edges),
                    ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
                        l, m+1, h, even_edges),
                    ForestedGraphComplex.ContractEdgesGO.generate_operator(
                        2*l-2+h, l, m+1, h, even_edges),
                    filterfunc
                ) for m in m_range])


create_forested_top_cohom_table(range(2,5), range(0,10), False, is_admissible)

# xxx = ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(
#                         2, 3, 2, False)
# aaa = xxx.get_matrix()
# print(aaa)
# vvv = xxx.domain
# bbb = vvv.get_basis()
# print(len(bbb))

# fil = basis_filter(V1)
# print(fil)

# ggg = V1.get_basis()[2]
# print(ggg)
# H = deleteUnmarkedEdges(ggg, V1)
# print("...")

# print(ggg.adjacency_matrix())
# print("...")
# print(H.adjacency_matrix())

# print(are_hairs_connected(ggg, V1))

# V1.display_basis_plots()

# cols = [i for i, ok in enumerate(fil) if ok]
# A = Dc1[:, cols]

# print("...")
# print(A)
