# Test the Lie bracket of ordinary graphs
import OrdinaryGraphComplex as OGC
import SpecialGraphs
from sage.all import *

even_edges = True

def apply_t(G):
    """Replace the last vertex of G by three vertices connected to it in a triangle.
    There must be exactly three edges incident to the last vertex. Each is connected to
    one of the new vertices.
    """
    n = G.order()
    last_vertex_neighbors = G.neighbors(n - 1)
    if len(last_vertex_neighbors) != 3:
        raise ValueError("The last vertex must have exactly three incident edges.") 
    G1 = Graph(G).copy()
    G1.delete_vertex(n - 1)
    new_vertices = [n - 1, n, n + 1]
    G1.add_vertices(new_vertices)
    for i in range(3):
        G1.add_edge(new_vertices[i], last_vertex_neighbors[i])
    G1.add_edges([(new_vertices[0], new_vertices[1]),
                  (new_vertices[1], new_vertices[2]),
                  (new_vertices[2], new_vertices[0])])
    return G1


# compute product of two tetrahedrons
T = SpecialGraphs.wheel_graph(3)
T2 = apply_t(T)
T3 = apply_t(T2)
T4 = apply_t(T3)
T5 = apply_t(T4)
D5 = SpecialGraphs.cube3d_graph() #SpecialGraphs.chain_of_diamonds_graph(2)
# print("Graph D5:", D5.graph6_string())
D5t = apply_t(D5)
# print("Graph T:", T.graph6_string())
# print("Graph T2:", T2.graph6_string())

# ll = OGC.LieBracket.lie_bracket_single(T, T, even_edges)
v = OGC.LieBracket.lie_bracket_single_vector(T, T4, even_edges)
v2 = OGC.LieBracket.lie_bracket_single_vector(T2, T3, even_edges)
v3 = OGC.LieBracket.lie_bracket_single_vector(T2, D5, even_edges)
v4 = OGC.LieBracket.lie_bracket_single_vector(T, D5t, even_edges)

V1 = OGC.OrdinaryGVS(T.order(), T.size() - T.order() + 1, even_edges)
# print(v)
# print(v2)
V = OGC.OrdinaryGVS(13, 9, even_edges)

# print nonzero entries
# for i in range(len(v)):
#     if v[i] != 0 or v2[i] != 0:
#         print(i, v[i], v2[i])

# for (G, sgn) in ll:
#     (g6, sgn2) = V.graph_to_canon_g6(G)
#     print(g6)
    # print(g6, sgn * sgn2)


print("Target dimension:",V.get_dimension())
# V.display_vector(v)

# check closedness of resulting vecfor under adjoint differential
op1 = OGC.ContractEdgesGO.generate_operator(13, 9, even_edges)
op2 = OGC.ContractEdgesGO.generate_operator(14, 9, even_edges)
print("Loading matrices ...")
D1 = op1.get_matrix()
D2 = op2.get_matrix()

print("Computing symmetry factor matrix...")
M = V.symmetry_factor_matrix()

vvv = vector(v)
vvv2 = vector(v2)
vvv3 = vector(v3)
vvv4 = vector(v4)
vs = M*vvv
vs2 = M*vvv2
vs3 = M*vvv3
vs4 = M*vvv4

print("Multiplying...")
w = D2.transpose() * vs
print("Sum of absolute values of w:", sum(abs(x) for x in w))
w2 = D2.transpose() * vs2
print("Sum of absolute values of w2:", sum(abs(x) for x in w2))
w3 = D2.transpose() * vs3
print("Sum of absolute values of w3:", sum(abs(x) for x in w3))
w4 = D2.transpose() * vs4
print("Sum of absolute values of w4:", sum(abs(x) for x in w4))


# # sanity check
# # ww = D2.transpose()  * vvv
# # print(ww)

# check whether not exact
r1 = D1.rank()
print("Rank D1:", r1)

D1a = D1.transpose().augment(vs)
# r2 = D1a.rank()
# print("Rank D1 augmented:", r2) # should be one larger

# r3 = D1.transpose().augment(vs).augment(vs2).rank()
# print("Rank D1 augmented twice:", r3) # should be two larger

r4 = D1.transpose().augment(vs).augment(vs2).augment(vs3).augment(vs4).rank()
print("Rank D1 augmented four times:", r4) # should be four larger
 
r5 = D1a.augment(vs2).augment(vs3).rank()
print("Rank D1 augmented thrice (v1,v2,v3):", r5) # should be three larger

r6 = D1a.augment(vs2).augment(vs4).rank()
print("Rank D1 augmented thrice (v1,v2,v4):", r6) # should be three larger

# D1a.change_ring(GF(1009))

# x = D1a.solve_right(vs)  # to check that it works
# print("Solution found:", x)