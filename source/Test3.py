import OrdinaryGraphComplex
import CHairyGraphComplex
from sage.all import *

vs = OrdinaryGraphComplex.OrdinaryGVS(6, 8, True)

print(vs.is_valid())
print(vs.exists_basis_file())
# print(vs.get_dimension())
# vs.build_basis()

gc = OrdinaryGraphComplex.OrdinaryGC(
    range(6, 7), range(8, 9), True, ["contract"])
print(gc.sum_vector_space.vs_list)
gc.build_basis()
print(vs.is_valid())
print(vs.exists_basis_file())
# print(vs.get_dimension())
