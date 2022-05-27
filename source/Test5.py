import unittest
import itertools
import logging
import Log
import TestGraphComplex
import OrdinaryGraphComplex
from sage.all import *
import StoreLoad


log_file = "Forested_Unittest.log"

# if __name__ == '__main__':
#     tt = OrdinaryGraphComplex.OrdinaryGC(
#         range(0, 15), range(0, 8), False, {'contract'})
#     tt.build_basis(ignore_existing_files=True, n_jobs=1, progress_bar=True)
#     tt.build_matrix(ignore_existing_files=True, n_jobs=1, progress_bar=True)


# transpose for testing
go = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(11, 9, False)
lst, shape = go._load_matrix_list()

a, b = shape
newshape = (b, a)

matrix_list = [(j, i, v) for (i, j, v) in lst]
matrix_list.sort()

# copied
data_type = "M"
(d, t) = newshape
stringList = []
stringList.append("%d %d %s" % (d, t, data_type))
for (i, j, v) in matrix_list:
    stringList.append("%d %d %d" % (i + 1, j + 1, v))
stringList.append("0 0 0")
StoreLoad.store_string_list(
    stringList, "/root/linbox-1.7.0/examples/testmatrix.txt")
