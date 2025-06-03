# check the cocycles of top degree in the odd graph complex
from sage.all import *
import OrdinaryGraphComplex
import KneisslerGC


k = 2
v = 8*k+2
l = 4*k+2

opD = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v, l, False)

KnV0 = KneisslerGC.KneisslerGVS(l, 0, False)
KnV2 = KneisslerGC.KneisslerGVS(l, 2, False)
KnV3 = KneisslerGC.KneisslerGVS(l, 3, False)

KnV0.build_basis()
KnV2.build_basis()
KnV3.build_basis()

b0 = KnV0.get_basis_g6()
b2 = KnV2.get_basis_g6()
b3 = KnV3.get_basis_g6()

bb = opD.domain.get_basis_g6()

D = opD.get_matrix()

x = D.right_kernel_matrix()
# print(x)
# print(x.nrows())
# print(x.ncols())
# print(bb)

print("barrel part")
for i in range(len(bb)):
    if x[0,i] != 0:
        g6 = bb[i]
        # print(f"g6: {g6}, x: {x[0,i]}")
        if g6 in b0:
            print(g6)
            
print("non-barrel part")
for i in range(len(bb)):
    if x[0,i] != 0:
        g6 = bb[i]
        # print(f"g6: {g6}, x: {x[0,i]}")
        if g6 in b3:
            print(g6)

print("rest")
for i in range(len(bb)):
    if x[0,i] != 0:
        g6 = bb[i]
        # print(f"g6: {g6}, x: {x[0,i]}")
        if not g6 in b3 and not g6 in b0:
            print(g6)