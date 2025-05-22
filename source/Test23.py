from sage.all import *
import OrdinaryGraphComplex

gg = Graph("Gqd`_[")
V = OrdinaryGraphComplex.OrdinaryGVS(8, 5, True)
D = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(8, 5, True)

l = D.operate_on(gg)

for (ggg, sgn) in l:
    print(ggg.graph6_string(), sgn)
