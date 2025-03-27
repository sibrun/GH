import unittest
import logging
import Log
import TestGraphComplex
import HairyGraphBiComplex


log_file = "HGBiC_Unittest.log"

d_range = range(13)
h_min_range = range(-12, -1)
even_edges = False
even_hairs = False


class GraphComplexTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [HairyGraphBiComplex.HairyCeEt1hBiGC(d_range, h_min_range, even_edges, even_hairs)]


class CohomologyTest(TestGraphComplex.CohomologyTest):
    def setUp(self):
        self.gc_list = [HairyGraphBiComplex.HairyCeEt1hBiGC(d_range, h_min_range, even_edges, even_hairs)]


class SquareZeroTest(TestGraphComplex.SquareZeroTest):
    def setUp(self):
        self.gc_list = [HairyGraphBiComplex.HairyCeEt1hBiGC(d_range, h_min_range, even_edges, even_hairs)]


class TestAcyclic(TestGraphComplex.TestAcyclic):
    def setUp(self):
        self.dif_list = HairyGraphBiComplex.HairyCeEt1hBiGC(d_range, h_min_range, even_edges, even_hairs).get_operator_list()


def suite():
    suite = unittest.TestSuite()
    suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(CohomologyTest('test_cohomology_functionality'))
    suite.addTest(SquareZeroTest('test_square_zero'))
    suite.addTest(TestAcyclic('test_acyclic'))
    return suite


if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite for hairy graph bi complex -----")
    runner = unittest.TextTestRunner()
    runner.run(suite())
