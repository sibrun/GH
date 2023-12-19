"""Provide an abstract class for a graph complex."""

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
        - sum_vector_space (GraphVectorSpace.SumVectorSpace): Vector space.
        - operator_collection_list (list(GraphOperator.OperatorMatrixCollection)): List of operators (differentials).
    """
    __metaclass__ = ABCMeta

    def __init__(self, sum_vector_space, operator_collection_list):
        """Initialize the underlying vector space and the list of operators."""
        self.sum_vector_space = sum_vector_space
        self.operator_collection_list = operator_collection_list

    @abstractmethod
    def __str__(self):
        """Unique description of the graph complex.

        :return: Unique description of the graph complex.
        :rtype: str
        """
        pass

    def get_vector_space(self):
        """Return the underlying vector space of the graph complex.

        :return: Vector space of the graph complex.
        :rtype: GraphVectorSpace.SumVectorSpace
        """
        return self.sum_vector_space

    def get_operator_list(self):
        """Return the list of operators (differential).

        :return: List of operators of the graph complex.
        :rtype: list(GraphOperator.OperatorMatrixCollection)
        """
        return self.operator_collection_list

    def plot_info(self):
        """Plot information about the underlying vector space an operator collections."""
        self.sum_vector_space.plot_info()
        for op_collection in self.operator_collection_list:
            op_collection.plot_info()

    def build_basis(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        """Build the basis of the underlying vector space.

        :param ignore_existing_files: Option to ignore existing basis files. Ignore existing files and
               rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
               (Default: False).
        :type ignore_existing_files: bool
        :param n_jobs: Option to compute the basis of the different sub vector spaces in parallel
               using n_jobs parallel processes (Default: 1).
        :type n_jobs: int
        :param progress_bar: Option to show a progress bar (Default: False). Only active if the basis of
               different sub vector spaces ar not built in parallel.
        :type progress_bar: bool
        :param info_tracker: Option to plot information about the sub vector spaces in a web page.
               Only active if basis not built in parallel processes (Default: False).
        :type info_tracker: bool
        """
        self.sum_vector_space.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                          progress_bar=progress_bar, info_tracker=info_tracker)

    def build_matrix(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        """Build the matrices for all operators of the graph complex.

        :param ignore_existing_files: Option to ignore existing matrix files. Ignore existing files and
               rebuild the operator matrices if True, otherwise skip rebuilding a matrix file if there exists already a
               matrix file (Default: False).
        :type ignore_existing_files: bool
        :param n_jobs: Option to build different matrices in parallel using n_jobs parallel processes (Default: 1).
        :type n_jobs: int
        :param progress_bar: Option to show a progress bar (Default: False).
               Only active if different matrices are not built in parallel.
        :type progress_bar: bool
        :param info_tracker: Option to plot information about the sub vector spaces in a web page (Default: False).
               Only active if different matrices are not built in parallel.
        :type info_tracker: bool
        """

        for op_collection in self.operator_collection_list:
            op_collection.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar,
                                       info_tracker=info_tracker)

    def square_zero_test(self):
        """For each differential of the graph complex, test whether it squares to zero.

        Return a dictionary to report the success of the test for each differential.

        :return: Dictionary (differential -> square zero test successful)
        :rtype: dict(differential -> bool)
        """
        test_dict = {}
        for dif in self.operator_collection_list:
            if isinstance(dif, GraphOperator.Differential):
                (triv_l, succ_l, inc_l, fail_l) = dif.square_zero_test()
                success = (succ_l > 0 and fail_l == 0)
                test_dict.update({dif: success})
        return test_dict

    def compute_rank(self, sage=None, linbox=None, rheinfall=None, ignore_existing_files=False, n_jobs=1,
                     info_tracker=False):
        """Compute the ranks of the operator matrices.

        :param sage: Use sage to compute the rank. Options: 'integer' (exact rank over the integers),
            'mod' (rank over a finite field, i.e. calculations modulo a prime number (Default: None).
        :type sage: str or list(str)
        :param linbox: Use linbox to compute the rank. Options: 'rational' (exact rank over the rational numbers),
               'mod' (rank over a finite field, i.e. calculations modulo a prime number (Default: None).
        :type linbox: str or list(str)
        :param rheinfall: Use rhainfall to compute the rank. Options: 'int64', 'mpq', 'mpz' (Default: None).
        :type rheinfall: str or list(str)
        :param sort_key: Sort the operator matrices to schedule the rank computation according to the sort key:
               'work_estimate', 'size', 'entries' (Default: 'size').
        :type sort_key: str
        :param ignore_existing_files: Option to ignore existing rank files. Ignore existing files and
               recompute the ranks if True, otherwise skip recomputing the rank if there exists already a
               rank file (Default: False).
        :type ignore_existing_files: bool
        :param n_jobs: Option to compute different ranks in parallel using
               n_jobs parallel processes (Default: 1).
        :type n_jobs: int
        :param info_tracker: Option to plot information about the operator matrices in a web page (Default: False).
               Only active if different ranks are not computed in parallel.
        :type info_tracker: bool

        .. seealso:: - http://www.linalg.org/
                    - https://github.com/linbox-team/linbox/blob/master/examples/rank.C
                    - https://github.com/riccardomurri/rheinfall/blob/master/src.c%2B%2B/examples/rank.cpp
        """
        if sage is None and linbox is None and rheinfall is None:
            raise ValueError("compute_rank: At least one rank computation method needs to be specified.")

        for op_collection in self.operator_collection_list:
            op_collection.compute_rank(sage=sage, linbox=linbox, rheinfall=rheinfall,
                                       ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                       info_tracker=info_tracker)

    def plot_cohomology_dim(self, to_html=False, to_csv=False, x_plots=2):
        """Plot the cohomology dimensions for each differential of the graph complex

        Plot the cohomology dimensions as plot and/or table associated with the differential.
        :param to_html: Option to generate a html file with a table of the cohomology dimensions (Default: False).
        :type to_html: bool
        :param to_csv: Option to generate a csv file with a table of the cohomology dimensions (Default: False).
        :type to_csv: bool
        :param x_plots: Number of plots on the x-axis (Default: 2).
        :type x_plots: int
        """
        for dif in self.operator_collection_list:
            if isinstance(dif, GraphOperator.Differential):
                dif.plot_cohomology_dim(
                    to_html=to_html, to_csv=to_csv, x_plots=x_plots)

    def export_cohomology_dim_for_web(self):
        """Writes the cohomology dimensions to a js file to be used on the website.
        """
        for dif in self.operator_collection_list:
            if isinstance(dif, GraphOperator.Differential):
                dif.export_cohomology_dim_for_web()

    def test_pairwise_anti_commutativity(self, commute=False):
        """Test pairwise anti-commutativity / commutativity of the op_collections of the graph complex.

        Return a dictionary to report the success of the test for each combination of operators.

        :param commute: If True test for commutativity, otherwise test for anti-commutativity
                (Default: False).
        :type commute: bool
        :return: Dictionary (pair of operators -> square zero test successful)
        :rtype: dict(tuple(GraphOperator.OperatorMatrixCollection, GraphOperator.OperatorMatrixCollection)  -> bool)
        """
        test_dict = {}
        for (op_collection1, op_collection2) in itertools.combinations(self.operator_collection_list, 2):
            (triv_l, succ_l, inc_l, fail_l) = self.test_anti_commutativity(
                op_collection1, op_collection2, commute=commute)
            success = (succ_l > 0 and fail_l == 0)
            test_dict.update({(op_collection1, op_collection2): success})
        return test_dict

    def test_anti_commutativity(self, op_collection1, op_collection2, commute=False, eps=Parameters.commute_test_eps):
        """Tests whether two operators (differentials) anti-commute.

        Tests whether the two operators (differentials) op_collection1 and op_collection2 anti-commute or commute.
        Searches all possible quadruples of operators and reports for how many of them the test was trivially successful
        (because at least two matrices are trivial), successful, inconclusive (because matrices are missing) or
        unsuccessful.

        :param op_collection1: First operator (differential).
        :type op_collection1: GraphOperator.OperatorMatrixCollection
        :param op_collection2: Second operator (differential).
        :type op_collection2:  raphOperator.OperatorMatrixCollection
        :param commute: If True test for commutativity, otherwise test for anti-commutativity (Default: False).
        :type commute: bool
        :param eps: Threshold for equivalence of matrices (Default: Parameters.commute_test_eps).
        :type eps: float
        :return: Tuple with the number of quadruples for which the (anti-)commutativity test was a
                 (trivial success, success, inconclusive, fail).
        :rtype: tuple(int, int, int, int)
        """
        case = 'anti-commutativity' if not commute else 'commutativity'
        print(' ')
        print('Test %s for %s and %s' %
              (case, str(op_collection1), str(op_collection2)))
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
                                res = self._test_anti_commutativity_for_quadruple(
                                    quadruple, commute=commute, eps=eps)
                                if res == 'triv':
                                    triv_l += 1
                                elif res == 'succ':
                                    print('success')
                                    succ_l += 1
                                elif res == 'inc':
                                    inc_l += 1
                                elif res == 'fail':
                                    print("%s test for %s: failed for the quadruples %s, %s, %s, %s" %
                                          (case, str(self), str(op1a), str(op1b), str(op2a), str(op2b)))
                                    logger.error("%s test for %s: failed for the quadruples %s, %s, %s, %s" %
                                                 (case, str(self), str(op1a), str(op1b), str(op2a), str(op2b)))
                                    fail.append(quadruple)
                                else:
                                    raise ValueError(
                                        'Undefined commutativity test result')

        fail_l = len(fail)
        print("trivial success: %d, success: %d, inconclusive: %d, failed: %d quadruples" %
              (triv_l, succ_l, inc_l, fail_l))
        logger.warning('Test %s for %s and %s' %
                    (case, str(op_collection1), str(op_collection2)))
        logger.warning("trivial success: %d, success: %d, inconclusive: %d, failed: %d quadruples" %
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
                else:
                    return 'fail'
            if (op1a.is_trivial() or op2b.is_trivial()) and (not (op2a.is_trivial() or op1b.is_trivial())):
                M1b = op1b.get_matrix()
                M2a = op2a.get_matrix()
                if Shared.matrix_norm(M1b * M2a) < eps:
                    return 'succ'
                else:
                    return 'fail'
            M1a = op1a.get_matrix()
            M2b = op2b.get_matrix()
            M1b = op1b.get_matrix()
            M2a = op2a.get_matrix()
            if Shared.matrix_norm(M2b * M1a + (-1 if commute else 1) * M1b * M2a) < eps:
                return 'succ'
            return 'fail'
        except StoreLoad.FileNotFoundError:
            return 'inc'
