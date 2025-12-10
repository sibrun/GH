# compute cohomology of forested 2 hair with spine graph complex

from sage.all import *
import ForestedGraphComplex
import networkx as nx

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
def filter_from_vs(V, filterfunc):
    return [filterfunc(g, V) for g in V.get_basis() ]

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

def cohom_formatted_sumvs(D1, D2, filterfunc):
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
    
    vs_fil = filter_from_sumvs(vs, filterfunc)
    out1_fil = filter_from_sumvs(D1.get_target(), filterfunc)

    # d = vs.get_dimension()
    d = sum( 1 for ok in vs_fil if ok)
    # print("dim", d)
    # print(vs_fil)

    r1 = 0
    r2 = 0

    vs2 = D2.get_domain()
    vs2_fil = filter_from_sumvs(vs2, filterfunc)
    out2_fil = filter_from_sumvs(D2.get_target(), filterfunc)
    d2 = sum( 1 for ok in vs2_fil if ok)


    
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
        # A2_fil = filter_cols(A2, vs2_fil)
        r2 = A2_fil.rank()
        # print("A2", A2.nrows(), A2.ncols(), r2)
        # print("A2_fil", A2_fil.nrows(), A2_fil.ncols(), r2)

    # exact or not?
    r_str = "" #f"({d}+{rc2}-{r1}-{r2},{d},{d2})"

    # iso string
    cohomdim = d-r1-r2

    return str(cohomdim) + r_str 

def cohom_formatted(D1, D2, filterfunc):
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
    
    vs_fil = filter_from_vs(vs, filterfunc)
    out1_fil = filter_from_vs(D1.get_target(), filterfunc)

    # d = vs.get_dimension()
    d = sum( 1 for ok in vs_fil if ok)
    # print("dim", d)
    # print(vs_fil)

    r1 = 0
    r2 = 0

    vs2 = D2.get_domain()
    vs2_fil = filter_from_vs(vs2, filterfunc)
    out2_fil = filter_from_vs(D2.get_target(), filterfunc)
    d2 = sum( 1 for ok in vs2_fil if ok)


    
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
        # A2_fil = filter_cols(A2, vs2_fil)
        r2 = A2_fil.rank()
        # print("A2", A2.nrows(), A2.ncols(), r2)
        # print("A2_fil", A2_fil.nrows(), A2_fil.ncols(), r2)

    # exact or not?
    r_str = "" #f"({d}+{rc2}-{r1}-{r2},{d},{d2})"

    # iso string
    cohomdim = d-r1-r2

    return str(cohomdim) + r_str 

def alwaystrue(G,V):
    return True

def create_forested_top_cohom_table(l_range, m_range, h, even_edges, filterfunc):
    s = ""
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
        
def create_forested_cohom_table(l_range, m_range, h, even_edges, filterfunc):
    s = ""
    for l in l_range:
        print("l =", l)
        print(  [cohom_formatted_sumvs(
                    ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
                        l, m, h, even_edges),
                    ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
                        l, m+1, h, even_edges),
                    filterfunc
                ) for m in m_range])
        
def create_forested_cohom_table_contract(l_range, m_range, h, even_edges, filterfunc):
    s = ""
    for l in l_range:
        print("l =", l)
        for m in m_range:
            print("max m =",m)
            maxv = 2*l-2+h
            print(  [cohom_formatted(
                    ForestedGraphComplex.ContractEdgesGO.generate_operator(
                        maxv-mm,l, mm, h, even_edges),
                    ForestedGraphComplex.ContractEdgesGO.generate_operator(
                        maxv-mm+1,l, mm+1, h, even_edges),
                    filterfunc
                ) for mm in range(m+1)])


def display_dimensions_forested(l_range, m_range, h, even_edges):
    s = ""
    for l in l_range:
        print("l =", l)
        for m in m_range:
            print("max m =",m)
            maxv = 2*l-2+h
            print(  [ForestedGraphComplex.ForestedGVS(
                        maxv-mm, l, m-mm, h, even_edges).get_dimension()
                    for mm in range(m+1)])
            

def filtered_dimension(V, filterfunc):
    return sum(1 for g in V.get_basis() if filterfunc(g, V))

def filtered_dimension_sumvs(Vsum, filterfunc):
    return sum( filtered_dimension(VV, filterfunc) for VV in Vsum.get_vs_list())

def display_dimensions_forested_filtered(l_range, m_range, h, even_edges, filterfunc):
    s = ""
    for l in l_range:
        print("l =", l)
        for m in m_range:
            print("max m =",m)
            maxv = 2*l-2+h
            print(  [filtered_dimension(ForestedGraphComplex.ForestedGVS(
                        maxv-mm, l, m-mm, h, even_edges), filterfunc)
                    for mm in range(m+1)])

def tadpole_and_paredge_free(G,V):
    n = V.n_vertices
    m = V.n_unmarked_edges
    for i in range(n, n+m):
        nh = G.neighbors(i)
        if len(nh) == 2:
            j = nh[0]
            k = nh[1]
            if G.has_edge(j, k):
                return False
        elif len(nh) == 1:
            return False
    return True

def convert_to_multigraph(G,V):
    GG = Graph(G, multiedges=True, loops=True)
    n = V.n_vertices
    m = V.n_unmarked_edges
    h = V.n_hairs
    for i in range(n, n+m):
        nh = G.neighbors(i)
        if len(nh) == 2:
            j = nh[0]
            k = nh[1]
            GG.add_edge(j, k)
            GG.delete_vertex(i)
        elif len(nh) == 1:
            j = nh[0]
            GG.add_edge(j, j)
            GG.delete_vertex(i)
        else:
            print("Unexpected number of neighbors:", len(nh))

    return GG
        
def is_3_edge_connected_fast_old(G):
    H = G.networkx_graph()
    return nx.edge_connectivity(H) >= 3

def is_3_edge_connected_fast(G):
    """
    Check whether an undirected Sage multigraph (with possible loops)
    is 3-edge-connected. Parallel edges are collapsed into weights.
    Loops are ignored because they do not affect cuts.
    Uses NetworkX Stoer–Wagner global min-cut.
    """

    # Build weighted simple graph for NetworkX
    H = nx.Graph()

    # Add all vertices
    for v in G.vertices():
        H.add_node(v)

    # Count multiplicities (ignore loops)
    weights = {}
    for u, v, _ in G.edges():
        if u == v:
            continue      # loops irrelevant
        if u > v:
            u, v = v, u
        weights[(u, v)] = weights.get((u, v), 0) + 1

    # Add weighted edges to NetworkX graph
    for (u, v), w in weights.items():
        H.add_edge(u, v, weight=w)

    # Too small → not 3-edge-connected
    if H.number_of_nodes() < 2:
        return False

    # Stoer–Wagner global min-cut
    cut_value, partition = nx.stoer_wagner(H, weight="weight")

    return cut_value >= 3

def is_edge_triconnected_x(G,V):
    MG = convert_to_multigraph(G,V)
    res = is_3_edge_connected_fast(MG) 
    # if not res and V.n_vertices >= 5:
    #     print("n", V.n_vertices)
    #     print("Not 3-edge-connected:")
    #     print(MG)
    #     print("Adjacency matrix:")
    #     print(MG.adjacency_matrix()     )
    #     G.plot().save("temp/graph_forest.png")
    #     MG.plot().save("temp/graph.png")
    #     exit(1)
    return res and tadpole_and_paredge_free(G,V)

def is_edge_triconnected(G,V):
    MG = convert_to_multigraph(G,V)
    res = is_3_edge_connected_fast(MG) 
    # if not res and V.n_vertices >= 5:
    #     print("n", V.n_vertices)
    #     print("Not 3-edge-connected:")
    #     print(MG)
    #     print("Adjacency matrix:")
    #     print(MG.adjacency_matrix()     )
    #     G.plot().save("temp/graph_forest.png")
    #     MG.plot().save("temp/graph.png")
    #     exit(1)
    return res 

def check_operator_sumvs(D, filterfunc):
    vs = D.get_domain()
    vs_fil = filter_from_sumvs(vs, filterfunc)
    out_fil = filter_from_sumvs(D.get_target(), filterfunc)

    A = D.get_matrix()
    # iteratre over non-zero entries
    for (i,j) in A.nonzero_positions():
        if (not out_fil[i]) and (vs_fil[j]):
            print("Error: matrix has non-zero entry mapping to filtered out basis element")
            print("i,j =", i, j)
            # print relevant graphs
            g_out = D.get_target().get_basis()[i]
            g_in = D.get_domain().get_basis()[j]
            g_out.plot().save("temp/graph_out.png")
            g_in.plot().save("temp/graph_in.png")
            
            mg_in = convert_to_multigraph(g_in, D.get_domain().get_vs_from_basis_index(j))
            mg_in.plot().save("temp/graph_in_multigraph.png")
            mg_out = convert_to_multigraph(g_out, D.get_target().get_vs_from_basis_index(i))
            mg_out.plot().save("temp/graph_out_multigraph.png")

            print(is_3_edge_connected_fast(mg_in), is_3_edge_connected_fast(mg_out))
            return False

def _to_weighted_simple_networkx(G):
    """
    Convert a Sage multigraph G (possibly with loops) into a weighted simple
    NetworkX Graph H. Parallel edges -> weight (multiplicity). Loops ignored.
    """
    H = nx.Graph()
    for v in G.vertices():
        H.add_node(v)

    weights = {}
    for u, v, _ in G.edges():
        if u == v:
            continue                # ignore loops (they don't affect cuts)
        # normalize order to treat (u,v) same as (v,u)
        if u > v:
            u, v = v, u
        weights[(u, v)] = weights.get((u, v), 0) + 1

    for (u, v), w in weights.items():
        H.add_edge(u, v, weight=w)

    return H


def is_essentially_4_edge_connected(G):
    """
    Return True iff the Sage multigraph G (possibly with loops) is
    essentially 4-edge-connected, i.e.
      - global edge-connectivity >= 3, and
      - there is no proper 3-edge-cut that isolates a single vertex.
    """
    H = _to_weighted_simple_networkx(G)

    # very small graphs cannot be essentially 4-edge-connected
    if H.number_of_nodes() < 3:
        return False

    # check global edge connectivity (Stoer-Wagner uses the 'weight' attribute)
    lambda_global, _ = nx.stoer_wagner(H, weight="weight")
    if lambda_global < 3:
        return False

    # Now search for a proper 3-cut that isolates a single vertex v.
    # For each v, test minimum cut between v and some t != v.
    # If any mincut has value 3 and the side containing v is exactly {v},
    # we've found a proper 3-cut.
    nodes = list(H.nodes())
    for i, v in enumerate(nodes):
        # early necessary check: if total weighted degree of v < 3 it cannot be essentially 4-edge-connected
        deg_v = sum(d.get('weight', 1) for _, d in H[v].items())
        if deg_v < 3:
            return False

        # For vertex v, try min-cuts to other vertices t.
        # We can stop as soon as we find one t that isolates v with cut value 3.
        for j, t in enumerate(nodes):
            if t == v:
                continue
            # use networkx.minimum_cut with capacity='weight'
            cut_value, (setA, setB) = nx.minimum_cut(H, v, t, capacity='weight')
            if cut_value == 3:
                # check whether the side containing v is the singleton {v}
                if (setA == {v}) or (setB == {v}):
                    return False

    # no proper 3-cut found
    return True

def is_edge_fourconnected(G,V):
    MG = convert_to_multigraph(G,V)
    res = is_essentially_4_edge_connected(MG) 
    return res

# create_forested_top_cohom_table(range(1,4), range(0,10), 0, False, is_admissible)
# create_forested_top_cohom_table(range(1,6), range(0,10), 0, False, alwaystrue)
# create_forested_top_cohom_table(range(1,6), range(0,10), 0, False, tadpole_and_paredge_free)
# create_forested_top_cohom_table(range(1,6), range(0,10), 0, False, is_edge_triconnected)
# create_forested_cohom_table(range(1,6), range(0,10), 0, False, alwaystrue)

# create_forested_cohom_table(range(1,6), range(0,10), 0, False, is_edge_triconnected)
# check_operator_sumvs(
#     ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(
#         4, 2, 0, False),
#     is_edge_triconnected
# )

# create_forested_top_cohom_table(range(1,6), range(0,10), 0, False, is_edge_fourconnected)
# create_forested_cohom_table(range(1,6), range(0,10), 0, False, is_edge_fourconnected)
# create_forested_top_cohom_table(range(1,6), range(0,10), 0, False, is_edge_fourconnected)


# display_dimensions_forested(range(5,6), range(0,10), 0, False)
# display_dimensions_forested_filtered(range(5,6), range(0,10), 0, False, is_edge_triconnected)
# display_dimensions_forested_filtered(range(5,6), range(0,10), 0, False, is_edge_fourconnected)

# create_forested_cohom_table_contract(range(1,6), range(0,6), 0, False, is_edge_triconnected)
# create_forested_cohom_table_contract(range(1,6), range(0,6), 0, False, is_edge_triconnected)
create_forested_cohom_table_contract(range(1,6), range(0,6), 0, False, is_edge_fourconnected)

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
