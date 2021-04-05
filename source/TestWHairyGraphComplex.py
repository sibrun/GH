import unittest
import itertools
import logging
import Log
import TestGraphComplex
import WHairyGraphComplex
from sage.all import *

log_file = "WHGC_Unittest.log"

# v_range = range(1, 8)
# l_range = range(0, 7)
# h_range = range(1, 8)
# edges_types = [True, False]
# hairs_types = [True, False]

# tt = WHairyGraphComplex.WHairyGraphVS(4,4,1,2)
# tt.build_basis()
# #tt.plot_all_graphs_to_file()
# tt.display_basis_plots()

#tt = WHairyGraphComplex.WHairyGraphVS(4,3,0,1)
#print(tt.is_valid())
#tt.build_basis(ignore_existing_files=True)

#tt.plot_all_graphs_to_file()
#tt.display_basis_plots()

# WGC = WHairyGraphComplex.WHairyGC(range(0,6), range(0,5), range(1,3), range(2,3) , ['contract'])

# WGC.build_basis(progress_bar=False, info_tracker=False, ignore_existing_files=True)
# WGC.build_matrix(progress_bar=False, info_tracker=False, ignore_existing_files=True)

# # # WGC.build_basis(progress_bar=False, info_tracker=False, ignore_existing_files=False)
# # # WGC.build_matrix(progress_bar=False, info_tracker=False, ignore_existing_files=False)

# WGC.square_zero_test()



# tt = WHairyGraphComplex.WHairyGraphVS(3,4,1,2)
# tt.build_basis(ignore_existing_files=True)
# tt.plot_all_graphs_to_file(skip_existing=False)
# tt.display_basis_plots()

tt = WHairyGraphComplex.ContractEdgesGO.generate_operator(5,4,1,2)
tu = WHairyGraphComplex.ContractEdgesGO.generate_operator(4,4,1,2)
D1 = tt.get_matrix()
D2 = tu.get_matrix()
C = D2*D1
print(D1) 
print(D2) 
print(C) 
# tt.domain.plot_all_graphs_to_file(skip_existing=False)
# tt.domain.display_basis_plots()
# tu.domain.plot_all_graphs_to_file(skip_existing=False)
# tu.domain.display_basis_plots()
# tu.target.plot_all_graphs_to_file(skip_existing=False)
# tu.target.display_basis_plots()

v = C[2]
print(v)
for i,j in enumerate(v):
    if j != 0:
        print(i)

# g = Graph("EZq?")
# res = tt.operate_on(g)
# print(res)
# for G,v in res:
#     g6, s = tt.target.graph_to_canon_g6(G)
#     print(g6, v*s)


# class BasisTest(TestGraphComplex.BasisTest):
#     def setUp(self):
#         self.vs_list = [HairyGraphComplex.HairyGraphVS(v, l, h, even_edges, even_hairs) for (v, l, h, even_edges, even_hairs)
#                         in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


# class OperatorTest(TestGraphComplex.OperatorTest):
#     def setUp(self):
#         self.op_list = [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, even_edges, even_hairs) for
#                         (v, l, h, even_edges, even_hairs) in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


# class GraphComplexTest(TestGraphComplex.GraphComplexTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]


# class CohomologyTest(TestGraphComplex.CohomologyTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]


# class SquareZeroTest(TestGraphComplex.SquareZeroTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]

# class AntiCommutativityTest(TestGraphComplex.AntiCommutativityTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, False, False, ['contract', 'et1h'])]

# def suite():
#     suite = unittest.TestSuite()
#     suite.addTest(BasisTest('test_basis_functionality'))
#     suite.addTest(BasisTest('test_compare_ref_basis'))
#     suite.addTest(OperatorTest('test_operator_functionality'))
#     suite.addTest(OperatorTest('test_compare_ref_op_matrix'))
#     suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
#     suite.addTest(CohomologyTest('test_cohomology_functionality'))
#     suite.addTest(SquareZeroTest('test_square_zero'))
#     suite.addTest(AntiCommutativityTest('test_anti_commutativity'))
#     return suite


# if __name__ == '__main__':
#     print("\n#######################################\n" + "----- Start test suite for hairy graph complex -----")
#     runner = unittest.TextTestRunner()
#     runner.run(suite())