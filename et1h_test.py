import HairyGraphComplex as HGC
from sage.all import *
import Shared as SH
reload(HGC)

o1 = HGC.EdgeToOneHairGO.generate_operator(8,5,4,True,False)
o2 = HGC.EdgeToOneHairGO.generate_operator(8,4,5,True,False)
o1.get_domain().build_basis()
o1.get_target().build_basis()
o2.get_target().build_basis()
o2.get_domain()==o1.get_target()
o1.build_matrix()
o2.build_matrix()
m1 = o1.get_matrix()
m2 = o2.get_matrix()
m = m2*m1
print(m)
print(SH.matrix_norm(m))
p=32189
m.change_ring(GF(p))
print(m.rank())

b1 = o1.get_domain().get_basis(g6=False)




