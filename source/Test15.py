import OrdinaryGraphComplex
import CHairyGraphComplex
import OrdinaryVariants

D = OrdinaryVariants.ContractEdgesGOTriconnected.generate_operator(12,10,False)

A = D.get_matrix()

print(A.size())