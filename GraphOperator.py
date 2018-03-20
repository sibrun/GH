from abc import ABCMeta, abstractmethod
import itertools
import pandas
import scipy.sparse as sparse
from scipy.sparse.linalg import aslinearoperator as aslinearoperator
from scipy.linalg.interpolative import estimate_rank as estimate_rank
from sage.all import *
import Log
import StoreLoad as SL
import ParallelProgress as PP
import Parameters
import Display

logger = Log.logger.getChild('graph_operator')


class OperatorMatrix(object):
    __metaclass__ = ABCMeta

    data_type = "M"

    def __init__(self, domain, target):
        if not self.is_match(domain, target):
            raise ValueError("Domain %s and target %s don't match to build the operator matrix %s"
                             % (str(domain), str(target), str(self)))
        self.domain = domain
        self.target = target

    def get_domain(self):
        return self.domain

    def get_target(self):
        return self.target

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def get_matrix_file_path(self):
        pass

    @abstractmethod
    def get_rank_file_path(self):
        pass

    def get_ref_matrix_file_path(self):
        pass

    def get_ref_rank_file_path(self):
        pass

    def is_valid(self):
        return True

    @abstractmethod
    def get_work_estimate(self):
        pass

    @abstractmethod
    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, n_jobs=1, progress_bar=True):
        pass

    @staticmethod
    @abstractmethod
    def is_match(domain, target):
        pass

    def exists_matrix_file(self):
        return os.path.isfile(self.get_matrix_file_path())

    def exists_rank_file(self):
        return os.path.isfile(self.get_rank_file_path())

    def delete_matrix_file(self):
        if os.path.isfile(self.get_matrix_file_path()):
            os.remove(self.get_matrix_file_path())

    def delete_rank_file(self):
        if os.path.isfile(self.get_rank_file_path()):
            os.remove(self.get_rank_file_path())

    def _store_matrix_list(self, matrixList, shape, data_type=data_type):
        (d, t) = shape
        stringList = []
        stringList.append("%d %d %s" % (d, t, data_type))
        for (domainIndex, targetIndex, data) in matrixList:
            stringList.append("%d %d %d" % (domainIndex + 1, targetIndex + 1, data))
        stringList.append("0 0 0")
        SL.store_string_list(stringList, self.get_matrix_file_path())

    def _load_matrix_list(self):
        if not self.exists_matrix_file():
            raise SL.FileNotFoundError("Cannot load matrix, No matrix file found for %s: " % str(self))
        stringList = SL.load_string_list(self.get_matrix_file_path())
        (d, t, data_type) = stringList.pop(0).split(" ")
        shape = (d, t) = (int(d), int(t))
        if d != self.domain.get_dimension() or t != self.target.get_dimension():
            raise ValueError("%s: Shape of matrix doesn't correspond to the vector space dimensions"
                             % str(self.get_matrix_file_path()))
        tail = map(int, stringList.pop().split(" "))
        if not tail == [0, 0, 0]:
            raise ValueError("%s: End line missing or matrix not correctly read from file"
                             % str(self.get_matrix_file_path()))
        matrixList = []
        for line in stringList:
            (i, j, v) = map(int, line.split(" "))
            if i < 1 or j < 1:
                raise ValueError("%s: Invalid matrix index: %d %d" % (str(self.get_matrix_file_path()), i, j))
            if i > d or j > t:
                raise ValueError("%s: Invalid matrix index outside matrix size:"
                                 " %d %d" % (str(self.get_matrix_file_path()), i, j))
            matrixList.append((i - 1, j - 1, v))
        return (matrixList, shape)

    def get_matrix_list(self):
        (matrixList, shape) = self._load_matrix_list()
        return matrixList

    def get_shifted_matrix_list(self, domain_start, target_start):
        matrixList = self.get_matrix_list()
        shiftedMatrixList = []
        for (i, j, v) in matrixList:
            shiftedMatrixList.append((i + domain_start, j + target_start, v))
        return shiftedMatrixList

    def get_matrix_shape(self):
        try:
            header = SL.load_line(self.get_matrix_file_path())
            (d, t, data_type) = header.split(" ")
            (d, t) = (int(d), int(t))
        except SL.FileNotFoundError:
            try:
                d = self.domain.get_dimension()
                t = self.target.get_dimension()
            except SL.FileNotFoundError:
                raise SL.FileNotFoundError("Matrix shape of %s unknown: "
                                           "Build matrix or domain and target basis first" % str(self))
        return (t, d)

    def get_matrix_entries(self):
        if not self.is_valid():
            return 0
        try:
            (matrixList, shape) = self._load_matrix_list()
            return len(matrixList)
        except SL.FileNotFoundError:
            raise SL.FileNotFoundError("Matrix entries unknown for %s: No matrix file" % str(self))

    def is_trivial(self):
        if not self.is_valid():
            return True
        (t, d) = self.get_matrix_shape()
        if t == 0 or d == 0:
            return True
        if self.get_matrix_entries() == 0:
            return True
        return False

    def get_matrix_transposed(self):
        if not self.is_valid():
            logger.warn("No matrix: %s is not valid" % str(self))
            (d ,t) = (self.domain.get_dimension(), self.target.get_dimension())
            entriesList = []
        else:
            (entriesList, shape) = self._load_matrix_list()
            (d, t) = shape
        M = matrix(ZZ, d, t, sparse=True)
        for (i, j, v) in entriesList:
            M.add_to_entry(i, j, v)
        return M

    def get_matrix(self):
        return self.get_matrix_transposed().transpose()

    def get_matrix_scipy_transposed(self):
        data = []
        row_ind = []
        col_ind = []
        (entriesList, shape) = self._load_matrix_list()
        for (r, c, d) in entriesList:
            row_ind.append(r)
            col_ind.append(c)
            data.append(d)
        M = sparse.csc_matrix((data, (row_ind, col_ind)), shape=shape, dtype='d')
        return M

    def compute_rank(self, mode='est', n_primes=2, primes=Parameters.primes, eps=Parameters.estimate_rank_eps):
        if (not self.is_valid()) or self.is_trivial():
            rank_dict = {'exact': 0}
        else:
            rank_dict = {}
            try:
                if mode in {'exact', 'mod_p', 'all'}:
                    M = self.get_matrix_transposed()
                if mode in {'exact', 'all'}:
                    rank_exact = M.rank()
                    rank_dict.update({'exact': rank_exact})
                if mode in {'mod_p', 'all'}:
                    n = min(n_primes, len(primes))
                    for p in primes[0:n]:
                        M.change_ring(GF(p))
                        rank_mod_p = M.rank()
                        info = 'mod_%d' % p
                        rank_dict.update({info: rank_mod_p})
                if mode in {'est', 'all'}:
                    rank_est = estimate_rank(aslinearoperator(self.get_matrix_scipy_transposed()), eps=eps)
                    rank_dict.update({'est': rank_est})
            except SL.FileNotFoundError:
                raise SL.FileNotFoundError("Cannot compute rank of %s: First build operator matrix" % str(self))
        return rank_dict

    def store_rank(self, mode='est', n_primes=2, ignore_existing_files=False, skip_if_no_matrix=True):
        if not ignore_existing_files and self.exists_rank_file():
            return
        try:
            rank_dict = self.compute_rank(mode=mode, n_primes=2)
        except SL.FileNotFoundError as error:
            if not skip_if_no_matrix:
                raise error
            else:
                logger.warn("Skip computing rank of %s, since matrix is not built" % str(self))
                return
        rank_list = [str(rank) + ' ' + mode for (mode, rank) in rank_dict.items()]
        SL.store_string_list(rank_list, self.get_rank_file_path())

    def get_matrix_rank(self):
        return 0
        if not self.is_valid():
            return 0
        try:
            (rank, info) = SL.load_string_list(self.get_rank_file_path())
        except SL.FileNotFoundError:
            raise SL.FileNotFoundError("Cannot load matrix rank, No rank file found for %s: " % str(self))
        return int(rank)

    def get_sort_value(self):
        try:
            entries = self.get_matrix_entries()
        except SL.FileNotFoundError:
            entries = Parameters.MAX_ENTRIES
        return entries

    def get_params_dict(self):
        return self.domain.get_params_dict()

    def get_info_dict(self):
        try:
            shape = self.get_matrix_shape()
        except SL.FileNotFoundError:
            shape = None
        try:
            entries = self.get_matrix_entries()
        except SL.FileNotFoundError:
            entries = None
        try:
            m_rank = self.get_matrix_rank()
        except SL.FileNotFoundError:
            m_rank = None
        return {'valid': self.is_valid(), 'shape': shape, 'entries': entries, 'rank': m_rank}


class Operator(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def operate_on(self, graph):
        """For G a graph returns a list of pairs (GG, x), such that (operator)(G) = sum x GG."""
        pass


class GraphOperator(Operator, OperatorMatrix):
    __metaclass__ = ABCMeta

    def __init__(self, domain, target):
        super(GraphOperator, self).__init__(domain, target)

    @classmethod
    def generate_op_matrix_list(cls, vector_space):
        vs_list = vector_space.get_vs_list()
        op_matrix_list = []
        for (domain, target) in itertools.product(vs_list, vs_list):
            if cls.is_match(domain, target):
                op_matrix_list.append(cls(domain, target))
        return op_matrix_list

    def __str__(self):
        return '<%s graph operator, domain: %s>' % (self.get_type(), str(self.domain))

    def is_valid(self):
        return self.domain.is_valid() and self.target.is_valid()

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, n_jobs=1, progress_bar=True):
        if not self.is_valid():
            return
        if not ignore_existing_files and self.exists_matrix_file():
            return
        try:
            domainBasis = self.domain.get_basis(g6=False)
        except SL.FileNotFoundError:
            if not skip_if_no_basis:
                raise SL.FileNotFoundError("Cannot build operator matrix of %s: "
                                           "First build basis of the domain %s" % (str(self), str(self.domain)))
            else:
                logger.warn("Skip building operator matrix of %s "
                             "since basis of the domain %s is not built" % (str(self), str(self.domain)))
                return
        try:
            targetBasis6 = self.target.get_basis(g6=True)
        except SL.FileNotFoundError:
            if not skip_if_no_basis:
                raise SL.FileNotFoundError("Cannot build operator matrix of %s: "
                                           "First build basis of the target %s" % (str(self), str(self.target)))
            else:
                logger.warn("Skip building operator matrix of %s "
                             "since basis of the target %s is not built" % (str(self), str(self.target)))
                return

        shape = (d, t) = (self.domain.get_dimension(), self.target.get_dimension())
        if d == 0 or t == 0:
            self._store_matrix_list([], shape)
            return

        lookup = {G6: j for (j, G6) in enumerate(targetBasis6)}

        desc = 'Build matrix: Domain: ' + str(self.get_params_dict())
        listOfLists = PP.parallel_common_progress(self._generate_matrix_list, list(enumerate(domainBasis)), lookup,
                                                  n_jobs=n_jobs, progress_bar=progress_bar, desc=desc)

        matrixList = list(itertools.chain.from_iterable(listOfLists))
        self._store_matrix_list(matrixList, shape)

    def _generate_matrix_list(self, domainBasisElement, lookup):
        (domainIndex, G) = domainBasisElement
        imageList = self.operate_on(G)
        canonImages = dict()
        for (GG, prefactor) in imageList:
            (GGcanon6, sgn1) = self.target.graph_to_canon_g6(GG)
            sgn0 = canonImages.get(GGcanon6)
            sgn0 = sgn0 if sgn0 is not None else 0
            canonImages.update({GGcanon6: (sgn0 + sgn1 * prefactor)})
        matrixList = []
        for (image, factor) in canonImages.items():
            if factor:
                targetIndex = lookup.get(image)
                if targetIndex is not None:
                    matrixList.append((domainIndex, targetIndex, factor))
        return matrixList


class BiOperatorMatrix(OperatorMatrix):
    def __init__(self, domain, target, op_collection1, op_collection2):
        super(BiOperatorMatrix, self).__init__(domain, target)
        self.op_collection1 = op_collection1
        self.op_collection2 = op_collection2

    def __str__(self):
        return '<Bi operator matrix on domain: %s, and %s, %s' \
               % (str(self.domain), str(self.op_collection1), str(self.op_collection2))

    def get_matrix_file_path(self):
        pass

    def get_rank_file_path(self):
        pass

    def get_work_estimate(self):
        return len(self.domain.get_vs_list())

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, n_jobs=1, progress_bar=True):
        if not ignore_existing_files and self.exists_matrix_file():
            return
        shape = (self.domain.get_dimensions(), self.target.get_dimensions())
        matrixList = []
        for op in self.op_collection1.get_op_list() + self.op_collection1.get_op_list():
            domain_start_idx = self.domain.get_start_idx(op.get_domain())
            target_start_idx = self.target.get_start_idx(op.get_target())
            subMatrixList = op.get_matrix_list()
            for (i, j, v) in subMatrixList:
                matrixList.append((i + domain_start_idx, j + target_start_idx, v))
        self._store_matrix_list(matrixList, shape)


class OperatorMatrixCollection(object):
    def __init__(self, vector_space, op_matrix_list):
        self.vector_space = vector_space
        self.op_matrix_list = op_matrix_list

    @abstractmethod
    def get_type(self):
        pass

    def __str__(self):
        return '<%s operator matrix collection>' % self.get_type()

    def get_op_list(self):
        return self.op_matrix_list

    def get_vector_space(self):
        return self.vector_space

    def sort(self, work_estimate=True):
        if work_estimate:
            self.op_matrix_list.sort(key=operator.methodcaller('get_work_estimate'))
        else:
            self.op_matrix_list.sort(key=operator.methodcaller('get_sort_value'))

    def build_matrix(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.plot_info()
        self.sort()
        for op in self.op_matrix_list:
            op.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar)
        self.plot_info()

    def compute_rank(self, mode='est', n_primes=2, ignore_existing_files=True, n_jobs=1):
        self.plot_info()
        self.sort(work_estimate=False)
        PP.parallel(self._compute_single_rank, self.op_matrix_list, n_jobs=n_jobs, mode=mode, n_primes=n_primes,
                    ignore_existing_files=ignore_existing_files)
        self.plot_info()

    def _compute_single_rank(self, op, mode='est', n_primes=2, ignore_existing_files=True):
        op.store_rank(mode=mode, n_primes=n_primes, ignore_existing_files=ignore_existing_files)

    def plot_info(self):
        opList = []
        for op in self.op_matrix_list:
            info_dict = op.get_info_dict()
            opList.append(op.get_params_dict().values() + [info_dict.get('valid'), info_dict.get('shape'),
                                                           info_dict.get('entries'), info_dict.get('rank')])
        opColumns = self.vector_space.get_params_range_dict().keys() + ['valid', 'shape', 'entries', 'rank']
        opTable = pandas.DataFrame(data=opList, columns=opColumns)
        opTable.sort_values(by=['valid', 'entries'], inplace=True, na_position='last')
        Display.display_pandas_df(opTable)


class Differential(OperatorMatrixCollection):
    __metaclass__ = ABCMeta

    def __init__(self, vector_space, op_matrix_list):
        super(Differential, self).__init__(vector_space, op_matrix_list)

    @staticmethod
    # Check whether opD.domain == opDD.target
    def is_match(opD, opDD):
        return opD.get_domain() == opDD.get_target()

    @staticmethod
    # Computes the cohomology dimension, i.e., dim(ker(D)/im(DD)) = dim(opD.domain) - rankD - rankDD
    def cohomology_dim(opD, opDD):
        try:
            dimV = opD.get_domain().get_dimension()
        except SL.FileNotFoundError:
            logger.warn("Cannot compute cohomology: First build basis for %s " % str(opD.get_domain()))
            return None
        if dimV == 0:
            return 0
        if opD.is_valid():
            try:
                rankD = opD.get_matrix_rank()
            except SL.FileNotFoundError:
                logger.warn("Cannot compute cohomology: Matrix rank not calculated for %s " % str(opD))
                return None
        else:
            rankD = 0
        if opDD.is_valid():
            try:
                rankDD = opDD.get_matrix_rank()
            except SL.FileNotFoundError:
                logger.warn("Cannot compute cohomology: Matrix rank not calculated for %s " % str(opDD))
                return None
        else:
            rankDD = 0
        cohomologyDim = dimV - rankD - rankDD
        if cohomologyDim < 0:
            raise ValueError("Negative cohomology dimension for %s" % str(opD.domain))
        return cohomologyDim

    # Computes the cohomology, i.e., ker(D)/im(DD)
    def get_general_cohomology_dim_dict(self):
        cohomology_dim = dict()
        for (opD, opDD) in itertools.product(self.op_matrix_list, self.op_matrix_list):
            if Differential.is_match(opD, opDD):
                dim = Differential.cohomology_dim(opD, opDD)
                cohomology_dim.update({opD.domain: dim})
        return cohomology_dim

    def get_cohomology_dim(self):
        cohomology_dim = self.get_general_cohomology_dim_dict()
        dim_dict = dict()
        for vs in self.vector_space.get_vs_list():
            dim_dict.update({vs.get_params_tuple(): cohomology_dim.get(vs)})
        return dim_dict

    def square_zero_test(self, eps=Parameters.square_zero_test_eps):
        succ = []  # holds pairs for which test was successful
        fail = []  # failed pairs
        triv = []  # pairs for which test trivially succeeded because at least one operator is the empty matrix
        inc = []  # pairs for which operator matrices are missing
        for (op1, op2) in itertools.product(self.op_matrix_list, self.op_matrix_list):
            if Differential.is_match(op2, op1):
                # A composable pair is found
                p = (op1, op2)
                if not (op1.is_valid() and op2.is_valid()):
                    triv.append(p)
                    continue
                try:
                    M1 = op1.get_matrix()
                    M2 = op2.get_matrix()
                except SL.FileNotFoundError:
                    logger.warn("Cannot test square zero: "
                                 "Operator matrix not built for %s or %s" % (str(op1), str(op2)))
                    inc.append(p)
                    continue
                if op1.is_trivial() or op2.is_trivial():
                    triv.append(p)
                    continue
                if sum(map(abs, (M2 * M1).list())) < eps:
                    succ.append(p)
                else:
                    fail.append(p)
        (triv_l, succ_l, inc_l, fail_l) = (len(triv), len(succ), len(inc), len(fail))
        logger.warn("Square zero test for %s: trivial success: "
                     "%d, success: %d, inconclusive: %d, failed: %d pairs" % (str(self), triv_l, succ_l, inc_l, fail_l))
        if inc_l:
            logger.warn("Square zero test for %s: inconclusive: %d paris" % (str(self), inc_l))
        for (op1, op2) in fail:
            logger.error("Square zero test for %s: failed for the pair %s, %s" % (str(self), str(op1), str(op2)))
        return (triv_l, succ_l, inc_l, fail_l)
