import unittest
import itertools
import TestGraphComplex
import WRHairyGraphComplex


log_file = "WRHGC_Unittest.log"

v_range = range(1, 5)
l_range = range(4)
h_range = range(1, 3)
w_range = range(1, 3)


class BasisTest(TestGraphComplex.BasisTest):
    def setUp(self):
        self.vs_list = [WRHairyGraphComplex.WRHairyGraphVS(v, l, h, w) for (v, l, h, w)
                        in itertools.product(v_range, l_range, h_range, w_range)]


class OperatorTest(TestGraphComplex.OperatorTest):
    def setUp(self):
        self.op_list = [WRHairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, w) for
                        (v, l, h, w) in itertools.product(v_range, l_range, h_range, w_range)]


class GraphComplexTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [ WRHairyGraphComplex.WRHairyGC(v_range, l_range, h_range, w_range, ['contract']) ]

class CohomologyTest(TestGraphComplex.CohomologyTest):
    def setUp(self):
        self.gc_list = [ WRHairyGraphComplex.WRHairyGC(v_range, l_range, h_range, w_range, ['contract']) ]


class SquareZeroTest(TestGraphComplex.SquareZeroTest):
    def setUp(self):
        self.gc_list = [ WRHairyGraphComplex.WRHairyGC(v_range, l_range, h_range, w_range, ['contract']) ]


def suite():
    suite = unittest.TestSuite()
    suite.addTest(BasisTest('test_basis_functionality'))
    #suite.addTest(BasisTest('test_compare_ref_basis'))
    suite.addTest(OperatorTest('test_operator_functionality'))
    #suite.addTest(OperatorTest('test_compare_ref_op_matrix'))
    suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(CohomologyTest('test_cohomology_functionality'))
    suite.addTest(SquareZeroTest('test_square_zero'))
    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite for wrhairy graph complex -----")
    runner = unittest.TextTestRunner()
    runner.run(suite())
