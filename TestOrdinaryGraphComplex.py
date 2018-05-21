from sage.all import *
import unittest
import itertools
import logging
import Log
import TestGraphComplex
import OrdinaryGraphComplex


log_file = "OGC_Unittest.log"

v_range = range(3, 10)
l_range = range(0, 10)
edges_types = [True, False]


class BasisTest(TestGraphComplex.BasisTest):
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


class OperatorTest(TestGraphComplex.OperatorTest):
    def setUp(self):
        self.op_list = [OrdinaryGraphComplex.ContractEdgesGO.generate_operator(v, l, even_edges) for (v, l, even_edges) in
                        itertools.product(v_range, l_range, edges_types)]


class GraphComplexTest(TestGraphComplex.GraphComplexTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, even_edges, ['contract']) for even_edges in edges_types]
        self.gc_list += [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, False, ['delete_e'])]


class CohomologyTest(TestGraphComplex.CohomologyTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, even_edges, ['contract']) for even_edges in edges_types]
        self.gc_list += [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, False, ['delete_e'])]


class SquareZeroTest(TestGraphComplex.SquareZeroTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, even_edges, ['contract']) for even_edges in edges_types]
        self.gc_list += [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, False, ['delete_e'])]


class AntiCommutativityTest(TestGraphComplex.AntiCommutativityTest):
    def setUp(self):
        self.gc_list = [OrdinaryGraphComplex.OrdinaryGC(v_range, l_range, False, ['contract', 'delete_e'])]


def suite():
    suite = unittest.TestSuite()
    suite.addTest(BasisTest('test_perm_sign'))
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
    logging.warn("\n#####################################\n" + "----- Start test suite for ordinary graph complex -----")
    runner = unittest.TextTestRunner()
    runner.run(suite())
