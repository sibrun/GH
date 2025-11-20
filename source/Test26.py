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
print(Dc1)


opc2 = ForestedGraphComplex.ContractEdgesGO.generate_operator(nv, g, m+1, 2, False)
opc2.build_matrix()
V2 = opc2.domain
Dc2 = opc2.get_matrix()
b2 = V2.get_basis()
print(Dc2)

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

fil = basis_filter(V1)
print(fil)

# ggg = V1.get_basis()[2]
# print(ggg)
# H = deleteUnmarkedEdges(ggg, V1)
# print("...")

# print(ggg.adjacency_matrix())
# print("...")
# print(H.adjacency_matrix())

# print(are_hairs_connected(ggg, V1))

# V1.display_basis_plots()
