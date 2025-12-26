# Test the Lie bracket of ordinary graphs
import OrdinaryGraphComplex as OGC
import SpecialGraphs
from sage.all import *

even_edges = True


# compute product of two tetrahedrons
T = SpecialGraphs.wheel_graph(3)

ll = OGC.LieBracket.lie_bracket_single(T, T, even_edges)
v = OGC.LieBracket.lie_bracket_single_vector(T, T, even_edges)

V1 = OGC.OrdinaryGVS(T.order(), T.size() - T.order() + 1, even_edges)
print(v)
V = OGC.OrdinaryGVS(7, 6, even_edges)

# for (G, sgn) in ll:
#     (g6, sgn2) = V.graph_to_canon_g6(G)
#     print(g6)
    # print(g6, sgn * sgn2)


print(V.get_dimension())
# V.display_vector(v)

# check closedness of resulting vecfor under adjoint differential
op1 = OGC.ContractEdgesGO.generate_operator(7, 6, even_edges)
op2 = OGC.ContractEdgesGO.generate_operator(8, 6, even_edges)

D1 = op1.get_matrix()
D2 = op2.get_matrix()

M = V.symmetry_factor_matrix()

vvv = vector(v)

w = D2.transpose() * M * vvv
print(w)
# sanity check
# ww = D2.transpose()  * vvv
# print(ww)

# check whether not exact
r1 = D1.rank()
print("Rank D1:", r1)
vs = M*vvv
r2 = D1.transpose().augment(vs).rank()
print("Rank D1 augmented:", r2) # should be one larger