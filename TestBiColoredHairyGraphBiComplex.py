import unittest
import logging
import Log
import TestGraphComplex
import BiColoredHairyGraphBiComplex


log_file = "HGC_Unittest.log"

d_range = range(4, 11)
h_a_min_range = range(-9, 1)
h_b_min_range = range(-9, 1)


class GraphComplexTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [BiColoredHairyGraphBiComplex.BiColoredHairyContractSplitBiGC(d_range, h_a_min_range,
                                                                                     h_b_min_range, False, True, True)]


class SquareZeroTest(TestGraphComplex.SquerZeroTest):
    def setUp(self):
        self.gc_list = [BiColoredHairyGraphBiComplex.BiColoredHairyContractSplitBiGC(d_range, h_a_min_range,
                                                                                     h_b_min_range, False, True, True)]

def suite():
    suite = unittest.TestSuite()
    #suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(SquareZeroTest('test_square_zero'))
    return suite


if __name__ == '__main__':
    Log.set_log_file(log_file)
    Log.set_log_level('warning')
    logging.warn("\n#######################################\n" + "----- Start test suite for bi colored hairy graph bi complex -----")
    runner = unittest.TextTestRunner()
    runner.run(suite())