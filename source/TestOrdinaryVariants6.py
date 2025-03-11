import OrdinaryVariants
import OrdinaryGraphComplex

l = 9
for v in range(20):
    V = OrdinaryVariants.OrdinaryGVSFull(v,l,False)
    W = OrdinaryVariants.OrdinaryGVSTriconnected(v,l,False)
    U = OrdinaryGraphComplex.OrdinaryGVS(v,l,False)
    print(l,v,":",V.get_dimension(),W.get_dimension(), U.get_dimension())