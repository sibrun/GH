import unittest
import itertools
import logging
import Log
import TestGraphComplex
import HairyGraphComplex


log_file = "HGC_Unittest.log"

v_range = range(1, 8)
l_range = range(0, 7)
h_range = range(1, 8)
edges_types = [True, False]
hairs_types = [True, False]


class BasisTest(TestGraphComplex.BasisTest):
    def setUp(self):
        self.vs_list = [HairyGraphComplex.HairyGraphVS(v, l, h, even_edges, even_hairs) for (v, l, h, even_edges, even_hairs)
                        in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


class OperatorTest(TestGraphComplex.OperatorTest):
    def setUp(self):
        self.op_list = [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, even_edges, even_hairs) for
                        (v, l, h, even_edges, even_hairs) in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


class GraphComplexTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
                        for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
        self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
                         even_edges in edges_types]


class CohomologyTest(TestGraphComplex.CohomologyTest):
    def setUp(self):
        self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
                        for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
        self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
                         even_edges in edges_types]


class SquareZeroTest(TestGraphComplex.SquareZeroTest):
    def setUp(self):
        self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
                        for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
        self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
                         even_edges in edges_types]

class AntiCommutativityTest(TestGraphComplex.AntiCommutativityTest):
    def setUp(self):
        self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, False, False, ['contract', 'et1h'])]

def suite():
    suite = unittest.TestSuite()
    suite.addTest(BasisTest('test_basis_functionality'))
    suite.addTest(BasisTest('test_compare_ref_basis'))
    suite.addTest(OperatorTest('test_operator_functionality'))
    suite.addTest(OperatorTest('test_compare_ref_op_matrix'))
    suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(CohomologyTest('test_cohomology_functionality'))
    suite.addTest(SquareZeroTest('test_square_zero'))
    suite.addTest(AntiCommutativityTest('test_anti_commutativity'))
    return suite


if __name__ == '__main__':
    Log.set_log_file(log_file)
    Log.set_log_level('warning')
    logging.warn("\n#######################################\n" + "----- Start test suite for hairy graph complex -----")
    runner = unittest.TextTestRunner()
    runner.run(suite())