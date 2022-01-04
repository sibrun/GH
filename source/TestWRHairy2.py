import unittest
import itertools
import logging
import Log
import TestGraphComplex
import WRHairyGraphComplex
import GraphOperator
import json
from sage.all import *


WGC = WRHairyGraphComplex.WRHairyGC(range(0, 14), range(
    1, 4), range(0, 2), range(1, 2), ['contract'])

# WGC = WRHairyGraphComplex.WRHairyGC(range(0, 14), range(
#     0, 1), range(4, 5), range(1, 2), ['contract'])

# WGC.build_basis(progress_bar=False, info_tracker=False,
#                 ignore_existing_files=True)
# WGC.build_matrix(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=True)

# WGC.build_basis(progress_bar=False, info_tracker=False,
#                 ignore_existing_files=True)
# WGC.build_matrix(progress_bar=False, info_tracker=False,
#                  ignore_existing_files=True)

# # WGC.square_zero_test()

# # WGC.compute_rank(ignore_existing_files=True, sage="mod")
# WGC.compute_rank(ignore_existing_files=True, sage="integer")
# WGC.plot_cohomology_dim(to_html=True)
# Euler char
# WGC.print_dim_and_eulerchar()
# WGC.print_cohomology_dim()

WGC.export_cohomology_dim_for_web()

# for dif in WGC.operator_collection_list:
#     dd = dif.get_cohomology_dim_dict()
#     print(dd)
#     s = "data_sources.push({ name : 'GH', description : 'GH computation', data: ["
#     for k, v in dd.items():
#         if v:
#             print(k, v)
#             s += str(list(k) + [v]) + ",\n"
#     s += "] });"

#     print(s)

#     with open('../web/wrhairy/GHdata.js', 'w') as outfile:
#         outfile.write(s)
# print(dif)
# if isinstance(dif, GraphOperator.Differential):
# dif.plot_cohomology_dim()

# WGC.plot_cohomology_dim()

# tt = WRHairyGraphComplex.WRHairyGraphVS(4, 3, 1, 1)
# ss = WRHairyGraphComplex.WRHairyGraphVS(5, 3, 1, 1)
# uu = WRHairyGraphComplex.WRHairyGraphVS(3, 3, 1, 1)

# dd1 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(4, 3, 1, 1)
# dd2 = WRHairyGraphComplex.ContractEdgesGO.generate_operator(5, 3, 1, 1)

# D1 = dd1.get_matrix()
# D2 = dd2.get_matrix()
# # tt.build_basis()
# # tt.plot_all_graphs_to_file(skip_existing=True)
# # tt.display_basis_plots()
# D1 = D1.dense_matrix().change_ring(QQ)
# D2 = D2.dense_matrix().change_ring(QQ)

# print(D1)
# print(D2)

# print("Ranks: ", D1.rank(), D2.rank())
# print(ss.get_dimension(), "-D2->", tt.get_dimension(), "-D1->", uu.get_dimension())


# nsp = D1.right_kernel()  # (basis='pivot')
# print(nsp)
# img = D2.transpose().image()
# print(img)

# HH = nsp/img
# print(HH)

# nsp2 = D2.kernel()
# img2 = D1.image()

# print(nsp2)
# print(img2)
