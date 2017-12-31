from sage.all import *
import OrdinaryGraphComplex as OGC
import Shared as SH
import RefData as REF
import scipy.sparse as sparse
import logging

reload(SH)
reload(REF)
reload(OGC)

log_path = SH.get_path_from_current('log', 'test.log')
SH.generate_path(log_path)
logging.basicConfig(filename=log_path, level=logging.WARN)

op = OGC.ContractGO.get_operator(9,6,True)
ref_op = REF.RefOperator(op)
ref_matrix = ref_op.get_matrix()
matrix = op.get_matrix()
logging.warn(matrix+ref_matrix)




