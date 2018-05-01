import unittest
import itertools
import logging
from sage.all import *
import StoreLoad as SL
import TestGraphComplex as TGC
import HairyGraphComplex as HGC
import Parameters


log_file = "HGC_Unittest.log"

v_range = range(3, 8)
l_range = range(3, 7)
h_range = range(3, 7)
edges_types = [True, False]
hairs_types = [True, False]


class HGCBasisTest(TGC.BasisTest):
    def setUp(self):
        self.vs_list = [HGC.HairyGraphSumVS(v, l, h, even_edges, even_hairs) for (v, l, h, even_edges, even_hairs)
                        in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]

    def tearDown(self):
        self.vs_list = None

    def test_perm_sign(self):
        logging.warn('----- Test permutation sign -----')
        logging.warn('NOT IMPLEMENTED')         #TODO: implement sign test for hairy graphs


class HGCOperatorTest(TGC.OperatorTest):
    def setUp(self):
        self.op_list = [HGC.ContractEdges.generate_operator(v, l, h, even_edges, even_hairs) for (v, l, h, even_edges, even_hairs)
                        in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]

    def tearDown(self):
        self.op_list = None


class HGCGraphComplexTest(TGC.GraphComplexTest):
    def setUp(self):
        self.gc_list = [HGC.HairyGC(v_range, l_range, h_range, even_edges, even_hairs) for (even_edges, even_hairs)
                        in itertools.product(edges_types, hairs_types)]

    def tearDown(self):
        self.gc_list = None


def suite():
    log_path = os.path.join(Parameters.log_dir, log_file)
    SL.generate_path(log_path)
    logging.basicConfig(filename=log_path, level=logging.WARN)
    logging.warn("\n#######################################\n" + "----- Start test suite for hairy graph complex -----")
    suite = unittest.TestSuite()
    suite.addTest(HGCBasisTest('test_perm_sign'))
    suite.addTest(HGCBasisTest('test_basis_functionality'))
    suite.addTest(HGCBasisTest('test_compare_ref_basis'))
    suite.addTest(HGCOperatorTest('test_operator_functionality'))
    suite.addTest(HGCOperatorTest('test_compare_ref_op_matrix'))
    suite.addTest(HGCGraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(HGCGraphComplexTest('test_compare_ref_cohomology'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())