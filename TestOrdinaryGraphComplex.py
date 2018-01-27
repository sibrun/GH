import unittest
import itertools
import logging
import OrdinaryGraphComplex as OGC
import Shared as SH
import RefData as REF

reload(OGC)
reload(SH)
reload(REF)

log_dir = "log"
log_file = "test.log"
test_file = "test.txt"
skip_existing_files = False

eps = 1.0e-6
n_jobs = 1


class OGCTestCase(unittest.TestCase):

    def test_perm_sign(self):
        logging.warn('----- Test permutation sign -----')

        vs_even = OGC.OrdinaryGVS(6, 5, even_edges=True)
        vs_odd = OGC.OrdinaryGVS(6, 5, even_edges=False)
        G = graphs.WheelGraph(5)
        p = [0, 2, 3, 4, 5, 1]
        self.assertTrue(vs_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        self.assertTrue(vs_even.perm_sign(G, p) == 1, 'incorrect permutation sign')
        p = [0, 1, 5, 4, 3, 2]
        #self.assertTrue(vs_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        #self.assertTrue(vs_even.perm_sign(G, p) == -1, 'incorrect permutation sign')
        p = [0, 1, 2, 4, 3, 5]
        self.assertTrue(vs_odd.perm_sign(G, p) == -1, 'incorrect permutation sign')
        self.assertTrue(vs_even.perm_sign(G, p) == 1, 'incorrect permutation sign')

    def test_basis_functionality(self):
        logging.warn('----- Test basis functionality -----')
        for even_edges in [True, False]:
            vs = OGC.OrdinaryGVS(9, 9, even_edges=even_edges)
            vs.delete_basis_file()
            self.assertFalse(vs.exists_basis_file(),'basis should have been deleted')
            self.assertRaises(SH.NotBuiltError, vs.get_basis)
            vs.build_basis()
            self.assertTrue(vs.exists_basis_file(),'basis should exist')
            basis_g6 = vs.get_basis()
            self.assertTrue(len(basis_g6) > 0,'%s: basis_g6 has size 0' % str(vs))
            self.assertIsInstance(basis_g6,list,'type of basis_g6 is not a list')
            for i in range(len(basis_g6)):
                self.assertIsInstance(basis_g6[i], basestring, "%s: type of basis_g6[%d] is not 'string'" % (str(vs), i))
            basis = vs.get_basis(g6=False)
            self.assertEqual(len(basis), len(basis_g6),'%s: basis_g6 and basis have not the same size' % str(vs))
            self.assertIsInstance(basis,list,'type of basis is not a list')
            for i in range(len(basis_g6)):
                self.assertEqual(type(basis[i]), type(Graph()), "%s: basis[%d] has not 'sage graph' type" % (str(vs), i))
                self.assertEqual(basis_g6[i], basis[i].graph6_string(),'%s: error: basis[%d].graph6_string is not equal to basis_g6[%d]' % (str(vs), i, i))

    def test_basis(self):
        logging.warn('----- Compare basis with reference -----')
        v_range = range(6,10)
        l_range = range(5,10)
        even_range = [True, False]
        vs_list = [OGC.OrdinaryGVS(v, l, even_edges) for (v, l, even_edges) in itertools.product(v_range, l_range, even_range)]
        for vs in vs_list:
            if not vs.valid:
                logging.info('%s: no basis test, since not valid' % str(vs))
                continue
            if not skip_existing_files:
                vs.delete_basis_file()
            vs.build_basis()
            logging.info("%s: %s" % (str(vs), vs.get_info()))
            basis_list = vs.get_basis(g6=True)
            dimension = vs.get_dimension()
            self.assertEqual(len(basis_list), dimension, '%s: basis dimension not consistant' % str(vs))
            basis_set=set(basis_list)
            self.assertTrue(len(basis_list) == len(basis_set), '%s: basis contains duplicates' % str(vs))
            ref_vs = REF.RefVectorSpace(vs)
            if not ref_vs.exists_basis_file():
                logging.warn('%s: no reference file for basis' % str(vs))
                continue
            ref_basis_set = set(ref_vs.get_basis_g6())
            self.assertSetEqual(basis_set, ref_basis_set, '%s: basis not equal reference basis' % str(vs))

    def test_operator_matrix(self):
        logging.warn('----- Test operator matrix -----')
        v_range = range(6,10)
        l_range = range(6,10)
        even_range = [True, False]

        vs_list = [OGC.OrdinaryGVS(v, l, even_edges) for (v, l, even_edges) in itertools.product(v_range, l_range, even_range)]
        op_list = [OGC.ContractGO.get_operator(v, l, even_edges) for (v, l, even_edges) in itertools.product(v_range, l_range, even_range)]

        for vs in vs_list:
            vs.build_basis()

        for op in op_list:
            if not op.valid:
                logging.info('%s: no operator test, since not valid' % str(op))
                continue
            if not op.exist_domain_target_files():
                logging.warn('%s: no operator test, domain or target not built' % str(op))
                continue
            if not skip_existing_files and op.exists_matrix_file():
                op.delete_matrix_file()
                self.assertFalse(op.exists_matrix_file(), '%s matrix file should have been deleted' % str(op))
            op.build_matrix(n_jobs=n_jobs)
            self.assertTrue(op.exists_matrix_file(), '%s matrix file should exist' % str(op))
            logging.info("%s: %s" % (str(op), op.get_info()))

            M = op.get_matrix()
            shape = (m, n) = (M.nrows(), M.ncols())
            rank = M.rank()
            self.assertEqual(op.get_matrix_shape(), shape, '%s: matrix shape not consistent' % str(op))
            self.assertEqual(shape, (op.target.get_dimension(), op.domain.get_dimension()), '%s: matrix shape not consistent with vector space dimensions' % str(op))

            if not skip_existing_files and op.exists_rank_file():
                op.delete_rank_file()
                self.assertFalse(op.exists_rank_file(), '%s rank file should have been deleted' % str(op))
            op.compute_rank()
            self.assertTrue(op.exists_rank_file(), '%s rank file should exist' % str(op))
            self.assertEqual(op.get_rank(), rank,'%s: inconsistent rank' % str(op))
            entries = op.get_matrix_entries()
            self.assertEqual(op.is_trivial(), m == 0 or n == 0 or entries == 0,'%s: triviality check wrong' % str(op))

            ref_op = REF.RefOperator(op)
            if not ref_op.exists_matrix_file():
                logging.warn('%s: no reference file for operator matrix' % str(op))
                continue
            ref_M = ref_op.get_matrix_wrt_ref()
            ref_shape = (ref_M.nrows(), ref_M.ncols())
            ref_rank = ref_M.rank()
            self.assertEqual(ref_shape, shape, '%s: shape of matrix and reference matrix not equal' % str(op))
            if not ref_op.exists_rank_file():
                logging.warn('%s: no reference rank file' % str(op))
            else:
                self.assertEqual(ref_rank, ref_op.get_rank(),'%s: inconsistent reference rank' % str(op))
            self.assertEqual(rank, ref_rank,'%s: rank and reference rank not equal' % str(op))
            logging.info("%s: matrix rank: %d, ref matrix rank: %d" % (str(op), rank, ref_rank))
            ref_M_transformed = ref_op.get_matrix()
            self.assertEqual(ref_M_transformed.rank(), ref_rank,'%s: transformed reference matrix has not same rank' % str(op))
            if op.domain.even_edges:
                ref_M_transformed = -ref_M_transformed               #TODO: sign error in transformation matrix for even edges
            self.assertTrue(M == ref_M_transformed, '%s: matrix and transformed reference matrix not equal' % str(op))

    def test_graph_complex(self):
        logging.warn('----- Test graph complex -----')
        v_range = range(6,9)
        l_range = range(6,9)
        even_range = [True]

        for even_edges in even_range:
            ogc = OGC.OrdinaryGC(v_range, l_range, even_edges)

            #ogc.build_basis(skip_existing_files=skip_existing_files)
            #ogc.build_operator_matrix(skip_existing_files=skip_existing_files, n_jobs=n_jobs)
            #(triv_l, succ_l, inc_l, fail_l) = ogc.square_zero_test(eps)
            #self.assertTrue(fail_l == 0, "%s: square zero test failed for %d pairs" % (str(ogc),fail_l))
            #ogc.compute_ranks(skip_existing_files=skip_existing_files)
            ogc.compute_cohomology()
            ogc.store_member_info()
            ogc.store_cohomology_dim()
            ogc.plot_cohomology_dim()


def suite():
    log_path = SH.get_path_from_current(log_dir, log_file)
    SH.generate_path(log_path)
    logging.basicConfig(filename=log_path, level=logging.WARN)
    logging.warn("###################################\n" + "----- Start test suite -----")
    suite = unittest.TestSuite()
    #suite.addTest(OGCTestCase('test_perm_sign'))
    #suite.addTest(OGCTestCase('test_basis_functionality'))
    #suite.addTest(OGCTestCase('test_basis'))
    suite.addTest(OGCTestCase('test_operator_matrix'))
    #suite.addTest(OGCTestCase('test_graph_complex'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())

