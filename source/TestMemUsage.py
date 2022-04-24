
import OrdinaryGraphComplex

import os, psutil
from time import perf_counter

process = psutil.Process(os.getpid())

op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(13,10,False)
vs = op.domain
# vs2=op.target
# lookup = {G6: j for (j, G6) in enumerate(vs2.get_basis_g6())}
# ppp=vs.get_partition()
# inimem = process.memory_info().rss
# tt = perf_counter()
# print(inimem)  # in bytes 
# for (i,G) in enumerate(vs.get_basis()):
#     if i==0:
#         inimem = process.memory_info().rss
#     # op.operate_on(G)
#     # op._generate_matrix_list2((i,G), lookup)
#     canonG, perm_dict = G.canonical_label(
#             partition=ppp, certificate=True, algorithm='sage')
#     # canonG = G.canonical_label(
#     #         partition=ppp, certificate=False)
#     if i%1000 == 0:
#         print(f"graph {i}, mem usage {process.memory_info().rss - inimem} time taken {perf_counter() - tt}")
#         tt = perf_counter()


print(process.memory_info().rss) 

lll = list(vs.get_basis())

print(process.memory_info().rss)  # in bytes 


