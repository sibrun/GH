import unittest
import itertools
import logging
from sage.all import *
import StoreLoad as SL
import TestGraphComplex as TGC
import OrdinaryGraphComplex as OGC
import Parameters


log_file = "OGC_Unittest.log"

v_range = range(4, 10)
l_range = range(4, 10)
edges_types = [True, False]


class OGCBasisTest(TGC.BasisTest):
    def setUp(self):
        self.vs_list = [OGC.OrdinaryGVS(v, l, even_edges) for (v, l, even_edges) in
                        itertools.product(v_range, l_range, edges_types)]

    def tearDown(self):
        self.vs_list = None

    def test_perm_sign(self):
        logging.warn('----- Test permutation sign -----')

        vs_even = OGC.OrdinaryGVS(6, 5, even_edges=True)
        vs_odd = OGC.OrdinaryGVS(6, 5, even_edges=False)
        G = graphs.WheelGraph(5)
        p = [0, 2, 3, 4, 5, 1]
        self.assertTrue(vs_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        self.assertTrue(vs_even.perm_sign(G, p) == 1, 'incorrect permutation sign')
        p = [0, 1, 5, 4, 3, 2]
        #self.assertTrue(vs_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')    #TODO: fix perm sign test
        #self.assertTrue(vs_even.perm_sign(G, p) == -1, 'incorrect permutation sign')
        p = [0, 1, 2, 4, 3, 5]
        self.assertTrue(vs_odd.perm_sign(G, p) == -1, 'incorrect permutation sign')
        self.assertTrue(vs_even.perm_sign(G, p) == 1, 'incorrect permutation sign')


class OGCOperatorTest(TGC.OperatorTest):
    def setUp(self):
        self.op_list = [OGC.ContractDOrdinary.get_operator(v, l, even_edges) for (v, l, even_edges) in
                        itertools.product(v_range, l_range, edges_types)]

    def tearDown(self):
        self.op_list = None


class OGCGraphComplexTest(TGC.GraphComplexTest):
    def setUp(self):
        self.gc_list = [OGC.OrdinaryGC(v_range, l_range, even_edges) for even_edges in edges_types]

    def tearDown(self):
        self.gc_list = None


def suite():
    log_path = os.path.join(Parameters.log_dir, log_file)
    SL.generate_path(log_path)
    logging.basicConfig(filename=log_path, level=logging.WARN)
    logging.warn("\n#####################################\n" + "----- Start test suite for ordinary graph complex -----")
    suite = unittest.TestSuite()
    suite.addTest(OGCBasisTest('test_perm_sign'))
    suite.addTest(OGCBasisTest('test_basis_functionality'))
    suite.addTest(OGCBasisTest('test_compare_ref_basis'))
    suite.addTest(OGCOperatorTest('test_operator_functionality'))
    suite.addTest(OGCOperatorTest('test_compare_ref_op_matrix'))
    suite.addTest(OGCGraphComplexTest('test_graph_complex_functionality'))
    suite.addTest(OGCGraphComplexTest('test_compare_ref_cohomology'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
