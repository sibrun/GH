from sage.all import *
import OrdinaryGraphComplex as OGC
import Shared as SH
import RefData as REF
import scipy.sparse as sparse
import logging

reload(SH)
reload(REF)
reload(OGC)

m = matrix(ZZ, 4, 3, sparse=True)
m.add_to_entry(1, 2, -1)





