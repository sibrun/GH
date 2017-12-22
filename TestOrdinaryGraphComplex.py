import unittest
import itertools
import logging
import scipy.sparse as sparse
import scipy.linalg.interpolative as interpolative
from sage.all import *
import OrdinaryGraphComplex as OGC
import Shared as SH

reload(OGC)
reload(SH)

log_dir = "log"
log_file = "test.log"
test_file = "test.txt"
skip_existing_files = False

eps = 1.0e-6

class OGCTestCase(unittest.TestCase):

    def test_perm_sign(self):
        logging.info('----- Test permutation sign -----')

        ogc_even = OGC.OrdinaryGVS(6, 5, even_edges=True)
        ogc_odd = OGC.OrdinaryGVS(6, 5, even_edges=False)
        G = graphs.WheelGraph(5)
        p = [0, 2, 3, 4, 5, 1]
        self.assertTrue(ogc_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        self.assertTrue(ogc_even.perm_sign(G, p) == 1, 'incorrect permutation sign')
        p = [0, 1, 5, 4, 3, 2]
        #self.assertTrue(ogc_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        #self.assertTrue(ogc_even.perm_sign(G, p) == -1, 'incorrect permutation sign')
        p = [0, 1, 2, 4, 3, 5]
        self.assertTrue(ogc_odd.perm_sign(G, p) == -1, 'incorrect permutation sign')
        self.assertTrue(ogc_even.perm_sign(G, p) == 1, 'incorrect permutation sign')

    def test_basis_functionality(self):
        logging.info('----- Test basis functionality -----')
        for even_edges in [True, False]:
            vs = OGC.OrdinaryGVS(9, 9, even_edges=even_edges)
            vs.delete_file()
            self.assertFalse(vs.basis_built(),'basis should have been deleted')
            self.assertRaises(SH.NotBuiltError, vs.get_basis)
            vs.build_basis()
            self.assertTrue(vs.basis_built(),'basis should exist')
            basis_g6=vs.get_basis()
            self.assertTrue(len(basis_g6) > 0,'%s: basis_g6 has size 0' % str(vs))
            self.assertIsInstance(basis_g6,list,'type of basis_g6 is not a list')
            for i in range(len(basis_g6)):
                self.assertIsInstance(basis_g6[i], basestring, "%s: type of basis_g6[%d] is not 'string'" % (str(vs), i))
            basis=vs.get_basis(g6=False)
            self.assertEqual(len(basis), len(basis_g6),'%s: basis_g6 and basis have not the same size' % str(vs))
            self.assertIsInstance(basis,list,'type of basis is not a list')
            for i in range(len(basis_g6)):
                self.assertEqual(type(basis[i]), type(Graph()), "%s: basis[%d] has not 'sage graph' type" % (str(vs), i))
                self.assertEqual(basis_g6[i], basis[i].graph6_string(),'%s: error: basis[%d].graph6_string is not equal to basis_g6[%d]' % (str(vs), i, i))

    def test_basis(self):
        logging.info('----- Compare basis with reference -----')
        v_range = range(6,11)
        l_range = range(5,11)
        even_range = [True,False]
        for (v, l, even) in itertools.product(v_range, l_range, even_range):
            vs = OGC.OrdinaryGVS(v, l, even_edges=even, header_ref=False)
            if not vs.valid:
                logging.info('%s: no basis test, since not valid' % str(vs))
                continue
            if skip_existing_files:
                vs.delete_file()
            vs.build_basis()
            basis_list=vs.get_basis(g6=True)
            basis_set=set(basis_list)
            self.assertTrue(len(basis_list) == len(basis_set), '%s: basis contains duplicates' % str(vs))
            if not vs.ref_file_available():
                logging.warn('%s: no reference file for basis' % str(vs))
                continue
            ref_basis_set = set(vs.get_ref_basis_g6())
            self.assertSetEqual(basis_set, ref_basis_set, '%s: basis not equal reference basis' % str(vs))

    def test_operator_matrix(self):
        logging.info('----- Test operator matrix -----')
        v_range = range(7,10)
        l_range = range(6,10)
        even_range = [True, False]
        eps = 1.0e-6
        for even_edges in [True, False]:
            gc = OGC.OrdinaryGC(v_range, l_range, even_edges, skip_existing_files=skip_existing_files)
            gc.build_basis()
            gc.build_operator()
            for op in gc.op_list:
                if not op.valid:
                    logging.info('%s: no operator test, since not valid' % str(op))
                    continue
                if not op.matrix_built():
                    logging.warn('%s: no operator test, since matrix not built' % str(op))
                    continue
                matrix = op.get_matrix()
                if not op.ref_file_available():
                    logging.warn('%s: no reference file for operator matrix' % str(op))
                    continue
                matrix_ref = op.get_ref_matrix()
                self.assertEqual(matrix.shape, matrix_ref.shape, "%s: matrix and reference matrix don't have the same shape" % str(op))
                self.assertEqual(matrix.nnz, matrix_ref.nnz, "%s: matrix and reference matrix don't have the same number of explicite vlaues" % str(op))
                #self.assertTrue(abs(sparse.linalg.norm(matrix) - sparse.linalg.norm(matrix_ref)) < eps, "%s: matrix and reference matrix don't have the same norm" % str(op))
                #if matrix.nnz:
                   # rank = interpolative.estimate_rank(sparse.linalg.aslinearoperator(matrix), eps)
                    #ref_rank = interpolative.estimate_rank(sparse.linalg.aslinearoperator(matrix_ref), eps)
                    #self.assertEqual(rank, ref_rank,"%s: estimated ranks for matrix and reference matrix not equal" % str(op))
                    #print(str(op) + ": rank of matrix: " + str(rank) + ", rank of reference matrix: " + str(ref_rank))

    def test_operator_square_zero(self):
        logging.info('----- Test operator matrix -----')
        v_range = range(7,11)
        l_range = range(6,11)
        for even_edges in [True, False]:
            gc = OGC.OrdinaryGC(v_range, l_range, even_edges, skip_existing_files=skip_existing_files)
            gc.build_basis()
            gc.build_operator()
            gc.store_member_list()
            #(triv_l, succ_l, inc_l, fail_l) = gc.square_zero_test(eps)
            #self.assertTrue(fail_l==0,'square zero test failed')


def suite():
    log_path = SH.get_path_from_current(log_dir, log_file)
    SH.generate_path(log_path)
    logging.basicConfig(filename=log_path, level=logging.WARNING)
    logging.info("----- Start test suite -----")
    suite = unittest.TestSuite()
    #suite.addTest(OGCTestCase('test_perm_sign'))
    #suite.addTest(OGCTestCase('test_basis_functionality'))
    #suite.addTest(OGCTestCase('test_basis'))
    suite.addTest(OGCTestCase('test_operator_matrix'))
    suite.addTest(OGCTestCase('test_operator_square_zero'))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())

