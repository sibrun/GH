from abc import ABCMeta, abstractmethod
import itertools
from tqdm import tqdm
import StoreLoad as SL
import Shared as SH
import Parameters
import Log

logger = Log.logger.getChild('graph_complex')


class GraphComplex(object):
    __metaclass__ = ABCMeta
    def __init__(self, sum_vector_space, operator_collection_list):
        self.sum_vector_space = sum_vector_space
        self.operator_collection_list = operator_collection_list

    @abstractmethod
    def __str__(self):
        pass

    def get_vector_space(self):
        return self.sum_vector_space

    def get_differential_list(self):
        return self.operator_collection_list

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False, info_tracker=False):
        self.sum_vector_space.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                          progress_bar=progress_bar, info_tracker=info_tracker)

    def build_matrix(self, ignore_existing_files=True, n_jobs=1, progress_bar=False, info_tracker=False):
        for dif in self.operator_collection_list:
            dif.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar,
                             info_tracker=info_tracker)

    def square_zero_test(self):
        for dif in self.operator_collection_list:
            dif.square_zero_test()

    def compute_rank(self, exact=False, n_primes=1, estimate=False, ignore_existing_files=True, n_jobs=1,
                     info_tracker=False):
        for dif in self.operator_collection_list:
            dif.compute_rank(exact=exact, n_primes=n_primes, estimate=estimate,
                             ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, info_tracker=info_tracker)

    def plot_cohomology_dim(self):
        for dif in self.operator_collection_list:
            dif.plot_cohomology_dim()

    def test_pairwise_commutativity(self, anti_commute=False):
        for (op_collection1, op_collection2) in itertools.combinations(self.operator_collection_list, 2):
            self.test_commutativity(op_collection1, op_collection2, anti_commute=anti_commute)

    def test_commutativity(self, op_collection1, op_collection2, anti_commute=False, eps=Parameters.commute_test_eps):
        print(' ')
        print('Test commutativity for %s and %s' % (str(op_collection1), str(op_collection2)))
        succ_l = 0  # number of quadruples for which test was successful
        fail = []  # failed quadruples
        triv_l = 0  # number of quadruples for which test trivially succeeded because at least one operator is the empty matrix
        inc_l = 0  # number of quadruples for which operator matrices are missing

        op_list1 = op_collection1.get_op_list()
        op_list2 = op_collection2.get_op_list()
        for (op1a, op2a) in itertools.product(op_list1, op_list2):
            if op1a.get_domain() == op2a.get_domain():
                for op1b in op_list1:
                    if op1b.get_domain() == op2a.get_target():
                        for op2b in op_list2:
                            if op2b.get_domain() == op1a.get_target() and op1b.get_target() == op2b.get_target():

                                quadruple = (op1a, op1b, op2a, op2b)
                                res = self._test_commutativity_for_quadruple(quadruple, anti_commute=anti_commute,
                                                                             eps=eps)
                                if res == 'triv':
                                    triv_l += 1
                                elif res == 'succ':
                                    print('success')
                                    logger.warn('success')
                                    succ_l += 1
                                elif res == 'inc':
                                    inc_l += 1
                                elif res == 'fail':
                                    print('fail')
                                    logger.error('fail')
                                    print("Commutativity test for %s: failed for the quadruples %s, %s, %s, %s" %
                                                 (str(self), str(op1a), str(op1b), str(op2a), str(op2b)))
                                    logger.error("Commutativity test for %s: failed for the quadruples %s, %s, %s, %s" %
                                                 (str(self), str(op1a), str(op1b), str(op2a), str(op2b)))
                                    fail.append(quadruple)
                                else:
                                    raise ValueError('Undefined commutativity test result')

        fail_l = len(fail)
        print("trivial success: %d, success: %d, inconclusive: %d, failed: %d quadruples" %
              (triv_l, succ_l, inc_l, fail_l))
        logger.warn('Test commutativity for %s and %s' % (str(op_collection1), str(op_collection2)))
        logger.warn("trivial success: %d, success: %d, inconclusive: %d, failed: %d quadruples" %
                    (triv_l, succ_l, inc_l, fail_l))
        return (triv_l, succ_l, inc_l, fail_l)


    def _test_commutativity_for_quadruple(self, quadruple, anti_commute=False, eps=Parameters.commute_test_eps):
        (op1a, op1b, op2a, op2b) = quadruple
        if not ((op1a.is_valid() and op2b.is_valid()) or (op2a.is_valid() and op1b.is_valid())):
            return 'triv'

        if op1a.is_valid() and op2b.is_valid() and (not (op2a.is_valid() and op1b.is_valid())):
            try:
                if op1a.is_trivial() or op2b.is_trivial():
                    return 'triv'
                M1a = op1a.get_matrix()
                M2b = op2b.get_matrix()
            except SL.FileNotFoundError:
                return 'inc'
            if SH.matrix_norm(M2b * M1a) < eps:
                return 'succ'
            return 'fail'

        if (not (op1a.is_valid() and op2b.is_valid())) and op2a.is_valid() and op1b.is_valid():
            try:
                if op1b.is_trivial() or op2a.is_trivial():
                    return 'triv'
                M1b = op1b.get_matrix()
                M2a = op2a.get_matrix()
            except SL.FileNotFoundError:
                return 'inc'
            if SH.matrix_norm(M1b * M2a) < eps:
                return 'succ'
            return 'fail'

        try:
            if (op1a.is_trivial() or op2b.is_trivial()) and (op2a.is_trivial() or op1b.is_trivial()):
                return 'triv'
            if (not (op1a.is_trivial() or op2b.is_trivial())) and (op2a.is_trivial() or op1b.is_trivial()):
                M1a = op1a.get_matrix()
                M2b = op2b.get_matrix()
                if SH.matrix_norm(M2b * M1a) < eps:
                    return 'succ'
                else: return 'fail'
            if (op1a.is_trivial() or op2b.is_trivial()) and (not (op2a.is_trivial() or op1b.is_trivial())):
                M1b = op1b.get_matrix()
                M2a = op2a.get_matrix()
                if SH.matrix_norm(M1b * M2a) < eps:
                    return 'succ'
                else: return 'fail'
            M1a = op1a.get_matrix()
            M2b = op2b.get_matrix()
            M1b = op1b.get_matrix()
            M2a = op2a.get_matrix()
            if SH.matrix_norm(M2b * M1a + (1 if anti_commute else -1) * M1b * M2a) < eps:
                return 'succ'
            else:
                return 'fail'
        except SL.FileNotFoundError:
            return 'inc'
