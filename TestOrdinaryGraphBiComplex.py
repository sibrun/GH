import unittest
import logging
import Log
import TestGraphComplex
import OrdinaryGraphBiComplex


log_file = "OGBiC_Unittest.log"

d_range = range(0, 17)
even_edges = False


class GraphComplexTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphBiComplex.OrdinaryCeDeleBiGC(d_range, even_edges)]


class CohomologyTest(TestGraphComplex.CohomologyTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphBiComplex.OrdinaryCeDeleBiGC(d_range, even_edges)]


class SquareZeroTest(TestGraphComplex.SquareZeroTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphBiComplex.OrdinaryCeDeleBiGC(d_range, even_edges)]


def suite():
    suite = unittest.TestSuite()
    suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(CohomologyTest('test_cohomology_functionality'))
    suite.addTest(SquareZeroTest('test_square_zero'))
    return suite


if __name__ == '__main__':
    Log.set_log_file(log_file)
    Log.set_log_level('warning')
    logging.warn("\n#######################################\n" + "----- Start test suite for ordinary graph bi complex -----")
    runner = unittest.TextTestRunner()
    runner.run(suite())