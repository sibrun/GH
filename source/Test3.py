import OrdinaryGraphComplex
import CHairyGraphComplex
import HairyGraphComplex
import ForestedGraphComplex
from sage.all import *
import os
# vs = OrdinaryGraphComplex.OrdinaryGVS(6, 8, True)

# print(vs.is_valid())
# print(vs.exists_basis_file())
# # print(vs.get_dimension())
# # vs.build_basis()

# gc = OrdinaryGraphComplex.OrdinaryGC(
#     range(6, 7), range(8, 9), True, ["contract"])
# print(gc.sum_vector_space.vs_list)
# gc.build_basis()
# print(vs.is_valid())
# print(vs.exists_basis_file())
# # print(vs.get_dimension())

# def runme(cmd):
#     ret = os.system(cmd)
#     if ret != 0:
#         raise RuntimeError("Error return value")

# runme("genbgL 3 3")

# runme("genbgL 30 30 'tempwert.txt'")


# print("Blabla alive....")


# op = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(4, 0, 0, True)

# print(op.is_valid())

# print(op.domain.is_valid())
# print(op.target.is_valid())


vs = HairyGraphComplex.HairyGraphVS(12,5,8, True, True)

print(vs.get_basis_file_path())

print(vs.get_dimension())