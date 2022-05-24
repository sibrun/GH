import ForestedGraphComplex
import CHairyGraphComplex
from sage.all import *

import timeit

# vs = ForestedGraphComplex.PreForestedGVS(12,7,1,0)
for m in range(12):
    vs = ForestedGraphComplex.ForestedGVS(12,7,m,0, True)
    # for pvs in vs.get_required_prevs():
    #     pvs.build_basis()
    # vs.build_basis()
    print(m)
    print(vs.get_dimension())

    start = timeit.timeit()
    print( sum(1 for G in vs.get_basis() if G.is_biconnected() ) )
    end = timeit.timeit()
    print(end - start)

    start = timeit.timeit()
    print( sum(1 for G in vs.get_basis() if len(list(G.bridges()))==0 ) )
    end = timeit.timeit()
    print(end - start)

    # b2 = [G for G in b1 if G.is_biconnected() ]

    # b3 = [G for G in b1 if len(list(G.bridges()))==0 ]

    # print(len(b2), len(b3))