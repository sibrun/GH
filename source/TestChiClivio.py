import OrdinaryGraphComplex
from sage.all import *

# chi_2i


def chigraph(i):
    G = Graph(2*i)
    for j in range(2*i):
        G.add_edge(j, (j+1) % (2*i))
    for j in range(i):
        G.add_edge(j, (j+i) % (2*i))
    return G


i = 7
G = chigraph(i)
# G.show()

V = OrdinaryGraphComplex.OrdinaryGVS(2*i, i+1, True)
B = V.get_basis_g6()
s, c = V.graph_to_canon_g6(G)

# autom_list = G.automorphism_group().gens()
# print(V._has_odd_automorphisms(G, autom_list))

# for p in autom_list:
#     pd = p.dict()
#     pp = [pd[j] for j in range(G.order())]
#     if V.perm_sign(G, pp) == -1:
#         print(pp)

print(s, c)
print(B)

ind = B.index(s)

opD = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(2*i, i+1, True)
D = opD.get_matrix()
n = D.ncols()
m = D.nrows()
la = matrix(1, n)
la[0, ind] = 1

b = matrix(m+1, 1)
b[m, 0] = 1
DD = D.stack(la)
# print(D)
#  DD)
v = DD.solve_right(b)
print(v)

print(",,", (DD*v-b))
