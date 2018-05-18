"""Module providing an abstract class for a graph complex."""

__all__ = ['GraphComplex']

from abc import ABCMeta, abstractmethod
import itertools
import GraphOperator
import StoreLoad
import Shared
import Parameters
import Log

logger = Log.logger.getChild('graph_complex')


class GraphComplex(object):
    """Graph complex.

    Attributes:
        sum_vector_space (GraphVectorSpace.SumVectorSpace): Vector space.

        operator_collection_list (list(GraphOperator.OperatorMatrixCollection)): List of operators (differentials).
    """
    __metaclass__ = ABCMeta

    def __init__(self, sum_vector_space, operator_collection_list):
        """Initialize the underlying vector space and the list of operators."""
        self.sum_vector_space = sum_vector_space
        self.operator_collection_list = operator_collection_list

    @abstractmethod
    def __str__(self):
        """Unique description of the graph complex.

        :return: str: Unique description of the graph complex.
        """
        pass

    def get_vector_space(self):
        """Returns the underlying vector space of the graph complex.

        :return: GraphVectorSpace.SumVectorSpace: Vector space of the graph complex.
        """
        return self.sum_vector_space

    def get_operator_list(self):
        """Returns the list of operators (differential).

        :return: list(GraphOperator.OperatorMatrixCollection): List of operators of the graph complex.
        """
        return self.operator_collection_list

    def build_basis(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        """Builds the basis of the underlying vector space.

        :param ignore_existing_files: bool, optional: Option to ignore existing basis files. Ignore existing files and
                rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
                (Default: False).
        :param n_jobs: positive int, optional: Option to compute the basis of the different sub vector spaces in parallel
                using n_jobs parallel processes (Default: 1).
        :param progress_bar: bool, optional: Option to show a progress bar (Default: False). Only active if the basis of
                different sub vector spaces ar not built in parallel.
        :param info_tracker: bool, optional: Option to plot information about the sub vector spaces in a web page.
                Only active if basis not built in parallel processes (Default: False).
        """
        self.sum_vector_space.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                          progress_bar=progress_bar, info_tracker=info_tracker)

    def build_matrix(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        """Builds the matrices for all operators of the graph complex.

        :param ignore_existing_files: bool, optional: Option to ignore  existing matrix files. Ignore existing files and
            rebuild the operator matrices if True, otherwise skip rebuilding a matrix file if there exists already a
            matrix file (Default: False).
        :param n_jobs: positive int, optional: Option to build different matrices in parallel using
                n_jobs parallel processes (Default: 1).
        :param progress_bar: bool, optional: Option to show a progress bar (Default: False).
            Only active if different matrices are not built in parallel.
        :param info_tracker: bool, optional: Option to plot information about the sub vector spaces in a web page
            (Default: False). Only active if different matrices are not built in parallel.
        """

        for op_collection in self.operator_collection_list:
            op_collection.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar,
                             info_tracker=info_tracker)

    def square_zero_test(self):
        """For each differential of the graph complex, tests whether it squares to zero."""
        for dif in self.operator_collection_list:
            if isinstance(dif, GraphOperator.Differential):
                dif.square_zero_test()

    def compute_rank(self, exact=False, n_primes=1, estimate=False, ignore_existing_files=False, n_jobs=1,
                     info_tracker=False):
        """Computes the ranks for all operators of the graph complex.

        :param exact: bool, optional: Compute the exact rank (Default: False)
        :param n_primes: non-negative int, optional: Number of primes. Determine the rank over a finite field w.r.t.
            different prime numbers (Default: 1). If set to 0 rank is not computed modulo a prime number.
        :param estimate: bool, optional: If True estimate rank using interpolative mthods offered by the scipy package
            (Default: False).
        :param sort_key: Sort the operator matrices to schedule the rank computation according to the sort key:
            'work_estimate', 'size', 'entries' (Default: 'size').
        :param ignore_existing_files: bool, optional: Option to ignore existing rank files. Ignore existing files and
            recompute the ranks if True, otherwise skip recomputing the rank if there exists already a
            rank file (Default: False).
        :param n_jobs: positive int, optional: Option to compute different ranks in parallel using
                n_jobs parallel processes (Default: 1).
        :param info_tracker: bool, optional: Option to plot information about the operator matrices in a web page
            (Default: False). Only active if different ranks are not computed in parallel.
        """
        for op_collection in self.operator_collection_list:
            op_collection.compute_rank(exact=exact, n_primes=n_primes, estimate=estimate,
                                       ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                       info_tracker=info_tracker)

    def plot_cohomology_dim(self, to_html=False, to_csv=False, x_plots=2):
        """Plots the cohomology dimensions for each differential of the graph complex

        Plot the cohomology dimensions as plot and/or table associated with the differential.
        :param to_html: bool, optional: Option to generate a html file with a table of the cohomology dimensions
            (Dafault: False).
        :param to_csv: bool, optional: Option to generate a csv file with a table of the cohomology dimensions
            (default: False).
        :param x_plots, optional: positive int: Number of plots on the x-axis (Default: 2).
        """
        for dif in self.operator_collection_list:
            if isinstance(dif, GraphOperator.Differential):
                dif.plot_cohomology_dim(to_html=to_html, to_csv=to_csv, x_plots=x_plots)

    def test_pairwise_anti_commutativity(self, commute=False):
        """Test pairwise anti-commutativity / commutativity of the op_collections of the graph complex

        :param commute (bool, optional): If True test for commutativity, otherwise test for anti-commutativity
                (Default: False).
        """
        for (op_collection1, op_collection2) in itertools.combinations(self.operator_collection_list, 2):
            self.test_anti_commutativity(op_collection1, op_collection2, commute=commute)

    def test_anti_commutativity(self, op_collection1, op_collection2, commute=False, eps=Parameters.commute_test_eps):
        """Tests whether two operators (differentials) anti-commute.

        Tests whether the two operators (differentials) op_collection1 and op_collection2 anti-commute or commute.
        Searches all possible quadruples of operators and reports for how many of them the test was trivially successful
        (because at least two matrices are trivial), successful, inconclusive (because matrices are missing) or
        unsuccessful.

        :param op_collection1: GraphOperator.OperatorMatrixCollection: First operator (differential).
        :param op_collection2: raphOperator.OperatorMatrixCollection: Second operator (differential).
        :param commute: ool, optional: If True test for commutativity, otherwise test for anti-commutativity
                (Default: False).
        :param eps: positive float, optional: Threshold for equivalence of matrices (Default: Parameters.commute_test_eps).
        :return: tuple(int, int, int, int): Tuple with the number of quadruples for which the (anti-)commutativity test was a
            (trivial success, success, inconclusive, fail).
        """
        case = 'anti-commutativity' if not commute else 'commutativity'
        print(' ')
        print('Test %s for %s and %s' % (case, str(op_collection1), str(op_collection2)))
        succ_l = 0  # number of quadruples for which test was successful
        fail = []  # failed quadruples
        triv_l = 0  # number of quadruples for which test trivially succeeded because at least two operator are the empty matrix
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
                                res = self._test_anti_commutativity_for_quadruple(quadruple, commute=commute, eps=eps)
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
                                    print("%s test for %s: failed for the quadruples %s, %s, %s, %s" %
                                                 (case, str(self), str(op1a), str(op1b), str(op2a), str(op2b)))
                                    logger.error("%s test for %s: failed for the quadruples %s, %s, %s, %s" %
                                                 (case, str(self), str(op1a), str(op1b), str(op2a), str(op2b)))
                                    fail.append(quadruple)
                                else:
                                    raise ValueError('Undefined commutativity test result')

        fail_l = len(fail)
        print("trivial success: %d, success: %d, inconclusive: %d, failed: %d quadruples" %
              (triv_l, succ_l, inc_l, fail_l))
        logger.warn('Test %s for %s and %s' % (case, str(op_collection1), str(op_collection2)))
        logger.warn("trivial success: %d, success: %d, inconclusive: %d, failed: %d quadruples" %
                    (triv_l, succ_l, inc_l, fail_l))
        return (triv_l, succ_l, inc_l, fail_l)

    def _test_anti_commutativity_for_quadruple(self, quadruple, commute=False, eps=Parameters.commute_test_eps):
        (op1a, op1b, op2a, op2b) = quadruple
        if not ((op1a.is_valid() and op2b.is_valid()) or (op2a.is_valid() and op1b.is_valid())):
            return 'triv'

        if op1a.is_valid() and op2b.is_valid() and (not (op2a.is_valid() and op1b.is_valid())):
            try:
                if op1a.is_trivial() or op2b.is_trivial():
                    return 'triv'
                M1a = op1a.get_matrix()
                M2b = op2b.get_matrix()
            except StoreLoad.FileNotFoundError:
                return 'inc'
            if Shared.matrix_norm(M2b * M1a) < eps:
                return 'succ'
            return 'fail'

        if (not (op1a.is_valid() and op2b.is_valid())) and op2a.is_valid() and op1b.is_valid():
            try:
                if op2a.is_trivial() or op1b.is_trivial():
                    return 'triv'
                M1b = op1b.get_matrix()
                M2a = op2a.get_matrix()
            except StoreLoad.FileNotFoundError:
                return 'inc'
            if Shared.matrix_norm(M1b * M2a) < eps:
                return 'succ'
            return 'fail'

        try:
            if (op1a.is_trivial() or op2b.is_trivial()) and (op2a.is_trivial() or op1b.is_trivial()):
                return 'triv'
            if (not (op1a.is_trivial() or op2b.is_trivial())) and (op2a.is_trivial() or op1b.is_trivial()):
                M1a = op1a.get_matrix()
                M2b = op2b.get_matrix()
                if Shared.matrix_norm(M2b * M1a) < eps:
                    return 'succ'
                else: return 'fail'
            if (op1a.is_trivial() or op2b.is_trivial()) and (not (op2a.is_trivial() or op1b.is_trivial())):
                M1b = op1b.get_matrix()
                M2a = op2a.get_matrix()
                if Shared.matrix_norm(M1b * M2a) < eps:
                    return 'succ'
                else: return 'fail'
            M1a = op1a.get_matrix()
            M2b = op2b.get_matrix()
            M1b = op1b.get_matrix()
            M2a = op2a.get_matrix()
            if Shared.matrix_norm(M2b * M1a + (-1 if commute else 1) * M1b * M2a) < eps:
                return 'succ'
            return 'fail'
        except StoreLoad.FileNotFoundError:
            return 'inc'
