import unittest
import itertools
import logging
from sage.all import *
import Log
import TestGraphComplex
import OrdinaryGraphComplex


log_file = "OGC_Unittest.log"

v_range = range(4, 10)
l_range = range(4, 9)
edges_types = [True, False]


class OGCBasisTest(TestGraphComplex.BasisTest):
    def setUp(self):
        self.vs_list = [OrdinaryGraphComplex.OrdinaryGVS(v, l, even_edges) for (v, l, even_edges) in
                        itertools.product(v_range, l_range, edges_types)]

    def test_perm_sign(self):
        logging.warn('----- Test permutation sign -----')

        vs_even = OrdinaryGraphComplex.OrdinaryGVS(6, 5, even_edges=True)
        vs_odd = OrdinaryGraphComplex.OrdinaryGVS(6, 5, even_edges=False)
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


class OGCOperatorTest(TestGraphComplex.OperatorTest):
    def setUp(self):
        self.op_list = [OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v, l, even_edges) for (v, l, even_edges) in
                        itertools.product(v_range, l_range, edges_types)]


class OGCTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, False, ['contract', 'delete_e'])]

def suite():
    Log.set_log_file(log_file)
    Log.set_log_level('warning')
    logging.warn("\n#####################################\n" + "----- Start test suite for ordinary graph complex -----")
    suite = unittest.TestSuite()
    #suite.addTest(OGCBasisTest('test_perm_sign'))
    #suite.addTest(OGCBasisTest('test_basis_functionality'))
    #suite.addTest(OGCBasisTest('test_compare_ref_basis'))
    #suite.addTest(OGCOperatorTest('test_operator_functionality'))
    #suite.addTest(OGCOperatorTest('test_compare_ref_op_matrix'))
    suite.addTest(OGCTest('test_graph_complex_functionality'))
    suite.addTest(OGCTest('test_graph_complex'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
