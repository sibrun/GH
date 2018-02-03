from abc import ABCMeta, abstractmethod
import unittest
from sage.all import *
import logging
import StoreLoad as SL
import RefData as REF


class BasisTest(unittest.TestCase):
    __metaclass__ = ABCMeta

    @abstractmethod
    def setUp(self):
        pass

    @abstractmethod
    def tearDown(self):
        pass

    def test_basis_functionality(self):
        logging.warn('----- Test basis functionality -----')
        for vs in self.vs_list:
            if not vs.valid:
                logging.info('%s: no basis test, since not valid' % str(vs))
                continue
            vs.delete_basis_file()
            self.assertFalse(vs.exists_basis_file(), 'basis should have been deleted')
            self.assertRaises(SL.NotBuiltError, vs.get_basis)
            vs.build_basis(ignore_existing_file=True)
            self.assertTrue(vs.exists_basis_file(), 'basis should exist')
            logging.info("%s: %s" % (str(vs), vs.get_info()))
            basis_g6 = vs.get_basis()
            if len(basis_g6) == 0:
                logging.warn('len(basis_g6) == 0 for %s' % str(vs))
            self.assertIsInstance(basis_g6, list, 'type of basis_g6 is not a list')
            for G6 in basis_g6:
                self.assertIsInstance(G6, basestring, "%s: type of basis_g6 element is not string" % str(vs))
            self.assertEqual(vs.get_dimension(), len(basis_g6), 'dimension not consistent for %s' % str(vs))
            self.assertEqual(len(basis_g6), len(set(basis_g6)), '%s: basis contains duplicates' % str(vs))
            basis = vs.get_basis(g6=False)
            self.assertEqual(len(basis), len(basis_g6), '%s: basis_g6 and basis have not the same size' % str(vs))
            self.assertIsInstance(basis, list, 'type of basis is not a list')
            for G in basis:
                self.assertEqual(type(G), type(Graph()), "%s: basis element has not 'sage graph' type" % str(vs))
            basis_g6_2 = [G.graph6_string() for G in basis]
            self.assertListEqual(basis_g6, basis_g6_2,
                             '%s: error: G.graph6_string is not equal to G6 for G in basis and G6 in basis_g6' % str(vs))

    def test_basis(self):
        logging.warn('----- Compare basis with reference -----')
        for vs in self.vs_list:
            if not vs.valid:
                logging.info('%s: no basis test, since not valid' % str(vs))
                continue
            vs.build_basis(ignore_existing_file=True)
            logging.info("%s: %s" % (str(vs), vs.get_info()))
            basis_g6 = vs.get_basis(g6=True)
            ref_vs = REF.RefVectorSpace(vs)
            if not ref_vs.exists_basis_file():
                logging.warn('%s: no basis test since no reference file for basis' % str(vs))
                continue
            ref_basis_g6 = ref_vs.get_basis_g6()
            ref_basis_set = set(ref_basis_g6)
            self.assertEqual(len(ref_basis_g6), len(ref_basis_set), '%s: reference basis contains duplicates' % str(vs))
            self.assertSetEqual(set(basis_g6), ref_basis_set, '%s: basis not equal reference basis' % str(vs))


class OperatorTest(unittest.TestCase):
    __metaclass__ = ABCMeta

    @abstractmethod
    def setUp(self):
        pass

    @abstractmethod
    def tearDown(self):
        pass

    def test_operator_functionality(self):
        logging.warn('----- Test operator functionality -----')
        for op in self.op_list:
            if not op.valid:
                logging.info('%s: no operator test, since not valid' % str(op))
                continue
            if not op.exist_domain_target_files():
                logging.warn('%s: no operator test, domain or target not built' % str(op))
                continue
            op.delete_matrix_file()
            self.assertFalse(op.exists_matrix_file(), '%s matrix file should have been deleted' % str(op))
            op.build_matrix(ignore_existing_file=True, n_jobs=1)
            self.assertTrue(op.exists_matrix_file(), '%s matrix file should exist' % str(op))
            logging.info("%s: %s" % (str(op), op.get_info()))
            M = op.get_matrix()
            shape = (m, n) = (M.nrows(), M.ncols())
            rank = M.rank()
            self.assertEqual(op.get_matrix_shape(), shape, '%s: matrix shape not consistent' % str(op))
            self.assertEqual(shape, (op.target.get_dimension(), op.domain.get_dimension()),
                             '%s: matrix shape not consistent with vector space dimensions' % str(op))
            op.delete_rank_file()
            self.assertFalse(op.exists_rank_file(), '%s rank file should have been deleted' % str(op))
            op.compute_rank()
            self.assertTrue(op.exists_rank_file(), '%s rank file should exist' % str(op))
            self.assertEqual(op.get_rank(), rank, '%s: inconsistent rank' % str(op))
            entries = op.get_matrix_entries()
            self.assertEqual(op.is_trivial(), m == 0 or n == 0 or entries == 0,
                             '%s: triviality check wrong' % str(op))
            op.delete_matrix_file()
            op.build_matrix(ignore_existing_file=True, n_jobs=4)
            M_parallel = op.get_matrix()
            self.assertTrue(M == M_parallel, '%s matrix not equal if computed with parallel jobs' % str(op))

    def test_operator_matrix(self):
        logging.warn('----- Test operator functionality -----')
        for op in self.op_list:
            if not op.valid:
                logging.info('%s: no operator test, since not valid' % str(op))
                continue
            if not op.exist_domain_target_files():
                logging.warn('%s: no operator test, domain or target not built' % str(op))
                continue
            op.build_matrix(ignore_existing_file=True, n_jobs=1)
            op.compute_rank(ignore_existing_file=True)
            M = op.get_matrix()
            shape = op.get_matrix_shape()
            rank = op.get_rank()
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
                continue
            else:
                self.assertEqual(ref_rank, ref_op.get_rank(), '%s: inconsistent reference rank' % str(op))
                self.assertEqual(rank, ref_rank, '%s: rank and reference rank not equal' % str(op))
                logging.info("%s: matrix rank: %d, ref matrix rank: %d" % (str(op), rank, ref_rank))
            ref_M_transformed = ref_op.get_matrix()
            self.assertEqual(ref_M_transformed.rank(), ref_rank,
                             '%s: transformed reference matrix has not same rank' % str(op))
            if op.domain.even_edges:
                ref_M_transformed = -ref_M_transformed      # TODO: sign error in transformation matrix for even edges
            self.assertTrue(M == ref_M_transformed,
                             '%s: matrix and transformed reference matrix not equal' % str(op))


class GraphComplexTest(unittest.TestCase):
    __metaclass__ = ABCMeta

    eps = 1.0e-6

    @abstractmethod
    def setUp(self):
        pass

    @abstractmethod
    def tearDown(self):
        pass

    def test_graph_complex_functionality(self):
        logging.warn('----- Test graph complex -----')
        for gc in self.gc_list:
            gc.build(ignore_existing_files=True, n_jobs=1)
            (triv_l, succ_l, inc_l, fail_l) = gc.square_zero_test(GraphComplexTest.eps)
            if succ_l == 0:
                logging.warn('%s: no successful pairs in square zero test' % str(gc))
            self.assertTrue(fail_l == 0, "%s: square zero test failed for %d pairs" % (str(gc),fail_l))

            gc.compute_ranks(ignore_existing_files=True)
            gc.compute_cohomology_dim()
            self.assertTrue(gc.exists_cohomology_file(), "%s: cohomology file should exist" % str(self))
            gc.plot_cohomology_dim()