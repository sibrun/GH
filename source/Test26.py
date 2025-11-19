# compute cohomology of forested 2 hair with spine graph complex

from sage.all import *
import ForestedGraphComplex

g = 2
nv = 2*g
m = 3

V1 = ForestedGraphComplex.ForestedGVS(nv, g, m, 2, False)
V1.build_basis()
b1 = V1.get_basis()
opc1 = ForestedGraphComplex.ContractEdgesGO.generate_operator(nv, g, m, 2, False)
opuu = ForestedGraphComplex.UnmarkEdgesGO.generate_operator(nv, g, m, 2, False)

Dc1 = opc1.get_matrix()


opc2 = ForestedGraphComplex.ContractEdgesGO.generate_operator(nv, g, m+1, 2, False)
V2 = opc2.domain()
Dc2 = opc2.get_matrix()
b2 = V2.get_basis()

def is_admissible(G):
    GG = G.copy()
    

def basis_filter(basis):
    return [1 if is_admissible(g) else 0 for g in basis]


