import unittest
import itertools
import logging
import Log
import TestGraphComplex
import ForestedGraphComplex
from sage.all import *

# VS = ForestedGraphComplex.ForestedGVS(6, 4, 3, 0, True)
# VS = ForestedGraphComplex.ForestedDegSlice(4, 3, 0, True)
# VS.build_basis()
# VS.display_basis_plots()

# FTD = ForestedGraphComplex.ContractUnmarkTopD(
#     range(6), range(15), range(1), True)

# FTD.build_matrix()
# FTD.compute_rank(sage="integer")

# print(FTD.get_cohomology_dim_dict())

# FTD.plot_cohomology_dim()

for even_edges in [False]:
    for h in range(4,6):
        for l in range(0,8-h):
            for m in range(2):
                # FD = ForestedGraphComplex.ForestedDegSlice(l,m,h, even_edges)
                # # print(FD.is_valid())
                # if FD.is_valid():
                #     for rep_index in range(len(Partitions(h))):
                #         print(f"Building projector {even_edges},{h},{l},{m}, {rep_index}...")
                #         pop = FD.get_isotypical_projector(rep_index)
                #         pop.build_matrix()
                #         pop.compute_rank(linbox="rational")
                op = ForestedGraphComplex.ContractUnmarkBiOM.generate_operator(l,m,h,even_edges)
                for rep_index in range(len(Partitions(h))):
                    print(f"Building projector {even_edges},{h},{l},{m}, {rep_index}...")
                    rop = op.restrict_to_isotypical_component(rep_index)
                    rop.build_matrix()
                    rop.compute_rank(linbox="rational")

