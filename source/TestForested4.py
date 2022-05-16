import ForestedGraphComplex
import CHairyGraphComplex
from sage.all import *

vs = ForestedGraphComplex.ForestedGVS(12,7,0,0)
b1 = vs.get_basis()

b2 = [G for G in b1 if G.is_biconnected() ]

b3 = [G for G in b1 if len(list(G.bridges()))>0]

print(len(b1), len(b2), len(b3))