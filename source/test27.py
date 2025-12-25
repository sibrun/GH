from sage.all import *
import OrdinaryGraphComplex


op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(18,10,False)
# op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(10,6,False)

D = op.get_matrix()



# V = D.right_kernel()

# print(V)
# p = 1009  # 4-digit prime
p = 100003  # 6-digit prime
D_finite = D.change_ring(GF(p))
V_finite = D_finite.right_kernel()
print(V_finite)