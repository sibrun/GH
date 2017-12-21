import unittest
import itertools
import logging
import pickle
import scipy.sparse as sparse
from scipy.linalg.interpolative import estimate_rank as estimate_rank
from sage.all import *
import OrdinaryGraphComplex as OGC
import Shared as SH

reload(OGC)
reload(SH)

log_dir = "./log"
log_file = "test.log"

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
            ogvs = OGC.OrdinaryGVS(9, 9, even_edges=even_edges)
            ogvs.delete_file()
            self.assertFalse(ogvs.basis_built(),'basis should have been deleted')
            self.assertRaises(SH.NotBuiltError, ogvs.get_basis)
            ogvs.build_basis()
            self.assertTrue(ogvs.basis_built(),'basis should exist')
            basis_g6=ogvs.get_basis()
            self.assertTrue(len(basis_g6)>0,'basis_g6 has size 0 for the parameters: ' + ogvs.params_to_string())
            self.assertIsInstance(basis_g6,list,'type of basis_g6 is not a list')
            for i in range(len(basis_g6)):
                self.assertIsInstance(basis_g6[i], basestring, 'type of basis_g6[i] is not string for the parameters: ' + ogvs.params_to_string())
            basis=ogvs.get_basis(g6=False)
            self.assertEqual(len(basis), len(basis_g6),' basis_g6 and basis have not the same size for the parameters: ' + ogvs.params_to_string())
            self.assertIsInstance(basis,list,'type of basis is not a list')
            for i in range(len(basis_g6)):
                self.assertEqual(type(basis[i]), type(Graph()), 'type of basis[i] is not sage graph for the parameters: ' + ogvs.params_to_string())
                self.assertEqual(basis_g6[i], basis[i].graph6_string(),'error: basis[i].graph6_string is not equal to basis_g6[i] for the parameters: ' + ogvs.params_to_string())

    def test_basis(self):
        logging.info('----- Compare basis with reference -----')
        v_range = range(5,9)
        l_range = range(4,9)
        even_range = [True,False]
        for (v, l, even) in itertools.product(v_range, l_range, even_range):
            ogvs = OGC.OrdinaryGVS(v, l, even_edges=even, header_ref=False)
            ogvs.delete_file()
            ogvs.build_basis()
            basis_list=ogvs.get_basis(g6=True)
            basis_set=set(basis_list)
            self.assertTrue(len(basis_list) == len(basis_set), 'Basis contains duplicates for the parameters: ' + ogvs.params_to_string())
            if not ogvs.ref_file_available():
                logging.warn('No basis reference file for the parameters: ' + ogvs.params_to_string())
                continue
            ref_basis_set = set(ogvs.get_ref_basis_g6())
            intersection = basis_set.intersection(ref_basis_set)
            self.assertTrue(len(intersection) == len(basis_set) == len(ref_basis_set), 'Basis not equal reference basis for the parameters: ' + ogvs.params_to_string())

    def test_operator_matrix(self):
        logging.info('----- Test operator matrix functionality -----')
        v_range = range(7,9)
        l_range = range(7,9)
        even_range = [True, False]
        eps_rank = 1.0e-6
        print(eps_rank)

        gc = OGC.OrdinaryGC(v_range, l_range, even_range, delete_old=True)
        gc.build_basis()
        gc.build_operator()
        for op in gc.op_list:
            if not op.matrix_built():
                continue
            matrix = op.get_matrix()
            if not op.ref_file_available():
                logging.warn('No operator reference file for the parameters: ' + op.params_to_string())
                continue
            matrix_ref = op.get_ref_matrix()
            self.assertEqual(matrix.shape, matrix_ref.shape, "Matrix and reference matrix don't have the same shape for the parameters: " + op.params_to_string())
            self.assertEqual(matrix.nnz, matrix_ref.nnz, "Matrix and reference matrix don't have the same number of explicite vlaues for the parameters: " + op.params_to_string())
            if matrix.nnz:
                rank = estimate_rank(sparse.linalg.aslinearoperator(matrix), eps_rank)
                ref_rank = estimate_rank(sparse.linalg.aslinearoperator(matrix_ref), eps_rank)
                self.assertEqual(rank, ref_rank,"Estimated ranks for matrix and reference matrix not equal for the parameters: " + op.params_to_string())
                print("rank: " + str(rank) + ", ref_rank: " + str(ref_rank) + ", for the parameters: " + op.params_to_string())

def suite():
    logging.basicConfig(filename=SH.get_log_path(log_dir, log_file), level=logging.INFO)
    logging.info("----- Start test suite -----")
    suite = unittest.TestSuite()
    #suite.addTest(OGCTestCase('test_perm_sign'))
    #suite.addTest(OGCTestCase('test_basis_functionality'))
    #suite.addTest(OGCTestCase('test_basis'))
    suite.addTest(OGCTestCase('test_operator_matrix'))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
