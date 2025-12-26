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
D5 = SpecialGraphs.cube3d_graph() #SpecialGraphs.chain_of_diamonds_graph(2)
print("Graph D5:", D5.graph6_string())
D5t = apply_t(D5)
# print("Graph T:", T.graph6_string())
# print("Graph T2:", T2.graph6_string())

# ll = OGC.LieBracket.lie_bracket_single(T, T, even_edges)
v = OGC.LieBracket.lie_bracket_single_vector(T, T3, even_edges)
v2 = OGC.LieBracket.lie_bracket_single_vector(T2, T2, even_edges)
v3 = OGC.LieBracket.lie_bracket_single_vector(T, D5, even_edges)

V1 = OGC.OrdinaryGVS(T.order(), T.size() - T.order() + 1, even_edges)
# print(v)
# print(v2)
V = OGC.OrdinaryGVS(11, 8, even_edges)

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
op1 = OGC.ContractEdgesGO.generate_operator(11, 8, even_edges)
op2 = OGC.ContractEdgesGO.generate_operator(12, 8, even_edges)

D1 = op1.get_matrix()
D2 = op2.get_matrix()

M = V.symmetry_factor_matrix()

vvv = vector(v)
vvv2 = vector(v2)
vvv3 = vector(v3)
vs3 = M*vvv3
print("Sum of absolute values of v:", sum(abs(x) for x in vvv))
print("Sum of absolute values of v2:", sum(abs(x) for x in vvv2))
print("Sum of absolute values of v3:", sum(abs(x) for x in vvv3))


w = D2.transpose() * M * vvv
print("Sum of absolute values of w:", sum(abs(x) for x in w))
w2 = D2.transpose() * M * vvv2
print("Sum of absolute values of w2:", sum(abs(x) for x in w2))
w3 = D2.transpose() * M * vvv3
print("Sum of absolute values of w3:", sum(abs(x) for x in w3))

# # sanity check
# # ww = D2.transpose()  * vvv
# # print(ww)

# check whether not exact
r1 = D1.rank()
print("Rank D1:", r1)
vs = M*vvv
D1a = D1.transpose().augment(vs)
r2 = D1a.rank()
print("Rank D1 augmented:", r2) # should be one larger
vs2 = M*vvv2
r3 = D1.transpose().augment(vs).augment(vs2).rank()
print("Rank D1 augmented twice:", r3) # should be two larger

r4 = D1.transpose().augment(vs2).rank()
print("Rank D1 augmented once otherwise (v1):", r4) # should be two larger

r5 = D1a.augment(vs2).augment(vs3).rank()
print("Rank D1 augmented thrice (v1,v2,v3):", r5) # should be three larger

# D1a.change_ring(GF(1009))

# x = D1a.solve_right(vs)  # to check that it works
# print("Solution found:", x)

# check D5 is indeed rep of missing class 
op0 = OGC.ContractEdgesGO.generate_operator(8, 5, even_edges)
D0 = op0.get_matrix()
r0 = D0.rank()
print("Rank D0:", r0)
V0 = OGC.OrdinaryGVS(8, 5, even_edges)
M0 = V0.symmetry_factor_matrix()
vd5_0 = vector(V0.graph_list_to_vector([(D5, 1)]))
vt3_0 = vector(V0.graph_list_to_vector([(T3, 1)]))

r0b = D0.transpose().augment(vd5_0).augment(vt3_0).rank()
print("Rank D0 augmented with D5 and T3:", r0b) # should be two larger
print(vd5_0)
print(vt3_0)