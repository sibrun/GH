import unittest
import logging
import Log
import TestGraphComplex
import HairyGraphBiComplex


log_file = "HGBiC_Unittest.log"

d_range = range(0, 13)
h_min_range = range(-12, -1)


class GraphComplexTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [HairyGraphBiComplex.HairyCeEt1hBiGC(d_range, h_min_range, False, False)]


class SquareZeroTest(TestGraphComplex.SquerZeroTest):
    def setUp(self):
        self.gc_list = [HairyGraphBiComplex.HairyCeEt1hBiGC(d_range, h_min_range, False, False)]


class TestAcyclic(TestGraphComplex.TestAcyclic):
    def setUp(self):
        self.dif_list = HairyGraphBiComplex.HairyCeEt1hBiGC(d_range, h_min_range, False, False).get_operator_list()


def suite():
    suite = unittest.TestSuite()
    suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(SquareZeroTest('test_square_zero'))
    suite.addTest(TestAcyclic('test_acyclic'))
    return suite


if __name__ == '__main__':
    Log.set_log_file(log_file)
    Log.set_log_level('warning')
    logging.warn("\n#######################################\n" + "----- Start test suite for hairy graph bi complex -----")
    runner = unittest.TextTestRunner()
    runner.run(suite())