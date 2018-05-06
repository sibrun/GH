from abc import ABCMeta, abstractmethod
import itertools
from tqdm import tqdm
import collections
import scipy.sparse as sparse
from scipy.sparse.linalg import aslinearoperator as aslinearoperator
from scipy.linalg.interpolative import estimate_rank as estimate_rank
from sage.all import *
import Log
import StoreLoad as SL
import Parallel
import Shared as SH
import Parameters
import PlotCohomology
import DisplayInfo

logger = Log.logger.getChild('graph_operator')


class OperatorMatrixProperties(object):
    def __init__(self):
        self.valid = None
        self.shape = None
        self.entries = None
        self.rank = None
        self.rank_mod_p = None
        self.rank_est = None

    @classmethod
    def names(cls):
        return ['valid', 'shape', 'entries', 'rank', 'rank_mod_p', 'rank_estimate']

    @staticmethod
    def sort_variables():
        return ['valid', 'entries']

    def list(self):
        return [self.valid, self.shape, self.entries, self.rank, self.rank_mod_p, self.rank_est]


class OperatorMatrix(object):
    __metaclass__ = ABCMeta

    data_type = "M"

    def __init__(self, domain, target):
        if not self.is_match(domain, target):
            raise ValueError("Domain %s and target %s don't match to build the operator matrix %s"
                             % (str(domain), str(target), str(self)))
        self.domain = domain
        self.target = target
        self.properties = OperatorMatrixProperties()

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
        return self.domain.is_valid() and self.target.is_valid()

    @abstractmethod
    def get_work_estimate(self):
        pass

    @abstractmethod
    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
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
        for (i, j, v) in matrixList:
            stringList.append("%d %d %d" % (i + 1, j + 1, v))
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

    def get_matrix_shape_entries(self):
        try:
            (matrixList, shape) = self._load_matrix_list()
            (d, t) = shape
            return ((t, d), len(matrixList))
        except SL.FileNotFoundError:
            raise SL.FileNotFoundError("Matrix shape and entries unknown for %s: No matrix file" % str(self))

    def get_matrix_entries(self):
        if not self.is_valid():
            return 0
        (shape, entries) = self.get_matrix_shape_entries()
        return entries

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
            logger.warn("Zero matrix: %s is not valid" % str(self))
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

    def compute_rank(self, exact=False, n_primes=1, primes=Parameters.primes, estimate=False,
                     eps=Parameters.estimate_rank_eps, ignore_existing_files=False, skip_if_no_matrix=True):
        if not self.is_valid():
            return
        if not ignore_existing_files:
            if self.exists_rank_file():
                return
        elif self.exists_rank_file():
            self.delete_rank_file()
        print('Compute matrix rank: Domain: ' + str(self.domain.get_ordered_param_dict()))
        try:
            rank_dict = self._compute_rank(exact=exact, n_primes=n_primes, primes=primes, estimate=estimate, eps=eps)
        except SL.FileNotFoundError as error:
            if skip_if_no_matrix:
                logger.info("Skip computing rank of %s, since matrix is not built" % str(self))
                return
            else:
                raise error
        self._store_rank_dict(rank_dict)

    def _compute_rank(self, exact=False, n_primes=1, primes=Parameters.primes, estimate=False,
                      eps=Parameters.estimate_rank_eps):
        if self.is_trivial():
            rank_dict = {'exact': 0}
        else:
            rank_dict = {}
            try:
                if exact:
                    M = self.get_matrix_transposed()
                    rank_exact = M.rank()
                    rank_dict.update({'exact': rank_exact})
                if n_primes >= 1:
                    n = min(n_primes, len(primes))
                    for p in primes[0:n]:
                        M = self.get_matrix_transposed()
                        M.change_ring(GF(p))
                        rank_mod_p = M.rank()
                        info = 'mod_%d' % p
                        rank_dict.update({info: rank_mod_p})
                if estimate and min(self.get_matrix_shape()) >= Parameters.min_size_for_rank_estimate:
                    rank_est = estimate_rank(aslinearoperator(self.get_matrix_scipy_transposed()), eps=eps)
                    if rank_est != min(self.get_matrix_shape()):
                        rank_est -= 1
                    rank_dict.update({'estimate': rank_est})
            except SL.FileNotFoundError:
                raise SL.FileNotFoundError("Cannot compute rank of %s: First build operator matrix" % str(self))
        return rank_dict

    def _store_rank_dict(self, update_rank_dict):
        try:
            rank_dict = self._load_rank_dict()
        except SL.FileNotFoundError:
            rank_dict = dict()
        rank_dict.update(update_rank_dict)
        rank_list = [str(rank) + ' ' + mode for (mode, rank) in rank_dict.items()]
        SL.store_string_list(rank_list, self.get_rank_file_path())

    def _load_rank_dict(self):
        if not self.is_valid():
            return {'exact': 0}
        try:
            rank_list = SL.load_string_list(self.get_rank_file_path())
        except SL.FileNotFoundError:
            raise SL.FileNotFoundError("Cannot load matrix rank, No rank file found for %s: " % str(self))
        rank_dict = dict()
        for line in rank_list:
            (rank, mode) = line.split(" ")
            rank_dict.update({mode: int(rank)})
        return rank_dict

    def _get_ranks(self):
        if not self.is_valid():
            return (0, 0, 0)
        rank_dict = self._load_rank_dict()
        rank_exact = rank_dict.pop('exact', None)
        if rank_exact == 0:
            return (0, 0, 0)
        rank_est = rank_dict.pop('estimate', None)
        rank_mod_p = None
        ranks_mod_p = rank_dict.values()
        if len(ranks_mod_p) >= 1:
            if len(set(ranks_mod_p)) == 1:
                rank_mod_p = ranks_mod_p[0]
            else:
                logger.warn('Ranks modulo different primes not equal for %s' % str(self))
            if rank_est is not None and rank_mod_p != rank_est:
                logger.warn('Rank modulo a prime and estimated rank not equal for %s' % str(self))
        return (rank_exact, rank_mod_p, rank_est)

    def get_matrix_rank(self):
        if not self.is_valid():
            return 0
        (rank_exact, rank_mod_p, rank_est) = self._get_ranks()
        if rank_exact is not None:
            return rank_exact
        if rank_mod_p is not None:
            return rank_mod_p
        logger.warn('Estimated rank for %s' % str(self))
        return rank_est

    def get_sort_size(self):
        try:
            sort_size = min(self.get_matrix_shape())
        except SL.FileNotFoundError:
            sort_size = Parameters.max_sort_value
        return sort_size

    def get_sort_entries(self):
        try:
            sort_entries = self.get_matrix_entries()
        except SL.FileNotFoundError:
            sort_entries = Parameters.max_sort_value
        return sort_entries

    def update_properties(self):
        self.properties.valid = self.is_valid()
        try:
            self.properties.shape = self.get_matrix_shape()
        except SL.FileNotFoundError:
            pass
        try:
            self.properties.entries = self.get_matrix_entries()
        except SL.FileNotFoundError:
            pass
        try:
            (self.properties.rank, self.properties.rank_mod_p, self.properties.rank_est) = self._get_ranks()
        except SL.FileNotFoundError:
            pass

    def get_properties(self):
        return self.properties


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
        op_matrix_list = []
        for (domain, target) in itertools.permutations(vector_space.get_vs_list(), 2):
            if cls.is_match(domain, target):
                op_matrix_list.append(cls(domain, target))
        return op_matrix_list

    def __str__(self):
        return '<%s graph operator, domain: %s>' % (self.get_type(), str(self.domain))

    def operate_on_list(self, graph_sgn_list):
        g6_coordinates_dict = self.target.get_g6_coordinates_dict()
        imageDict = dict()
        for (G1, sgn1) in graph_sgn_list:
            G1_image = self.operate_on(G1)
            for (G2, sgn2) in G1_image:
                (G2_g6, sgn_p) = self.target.graph_to_canon_g6(G2)
                coord = g6_coordinates_dict.get(G2_g6)
                if coord is None:
                    continue
                v = imageDict.get(coord)
                new_v = sgn1 * sgn2 * sgn_p
                if v is not None:
                    new_v += v
                if new_v:
                    imageDict.update({coord: new_v})
                elif v is not None:
                    imageDict.pop(coord)
        return imageDict.items()

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        if not self.is_valid():
            return
        if (not ignore_existing_files) and self.exists_matrix_file():
            return
        try:
            domainBasis = self.domain.get_basis(g6=False)
        except SL.FileNotFoundError:
            if not skip_if_no_basis:
                raise SL.FileNotFoundError("Cannot build operator matrix of %s: "
                                           "First build basis of the domain %s" % (str(self), str(self.domain)))
            else:
                logger.info("Skip building operator matrix of %s "
                             "since basis of the domain %s is not built" % (str(self), str(self.domain)))
                return
        try:
            targetBasis6 = self.target.get_basis(g6=True)
        except SL.FileNotFoundError:
            if not skip_if_no_basis:
                raise SL.FileNotFoundError("Cannot build operator matrix of %s: "
                                           "First build basis of the target %s" % (str(self), str(self.target)))
            else:
                logger.info("Skip building operator matrix of %s "
                             "since basis of the target %s is not built" % (str(self), str(self.target)))
                return

        shape = (d, t) = (self.domain.get_dimension(), self.target.get_dimension())
        if d == 0 or t == 0:
            self._store_matrix_list([], shape)
            return

        lookup = {G6: j for (j, G6) in enumerate(targetBasis6)}
        desc = 'Build matrix of %s operator: Domain: %s' % (str(self.get_type()), str(self.domain.get_ordered_param_dict()))

        #listOfLists = Parallel.parallel_common_progress(self._generate_matrix_list, list(enumerate(domainBasis)), lookup,
                                                  #n_jobs=n_jobs, progress_bar=progress_bar, desc=desc)

        #if not progress_bar:
        print(desc)
        listOfLists = []
        #for domainBasisElement in tqdm(list(enumerate(domainBasis)), desc=desc, disable=(not progress_bar)):
        for domainBasisElement in list(enumerate(domainBasis)):
            listOfLists.append(self._generate_matrix_list(domainBasisElement, lookup))
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
    __metaclass__ = ABCMeta

    def __init__(self, domain, target, operator_cls1, operator_cls2):
        super(BiOperatorMatrix, self).__init__(domain, target)
        self.operator_cls1 = operator_cls1
        self.operator_cls2 = operator_cls2

    @classmethod
    def generate_op_matrix_list(cls, graded_sum_vs, operator_cls1, operator_cls2):
        graded_sum_vs_list = graded_sum_vs.get_vs_list()
        bi_op_matrix_list = []
        for (domain, target) in itertools.permutations(graded_sum_vs_list, 2):
            if cls.is_match(domain, target):
                bi_op_matrix_list.append(cls(domain, target, operator_cls1, operator_cls2))
        return bi_op_matrix_list

    def __str__(self):
        return '<Bi operator matrix on domain: %s>' % str(self.domain)

    def is_valid(self):
        return True

    def get_work_estimate(self):
        return self.domain.get_dimension() * self.target.get_dimension()

    def build_matrix(self, ignore_existing_files=False, skip_if_no_matrices=False, progress_bar=False, **kwargs):
        if (not ignore_existing_files) and self.exists_matrix_file():
            return
        print(' ')
        print('Build matrix of %s' % str(self))
        shape = (self.domain.get_dimension(), self.target.get_dimension())
        underlying_matrices = self._get_underlying_matrices()
        self._build_underlying_matrices(underlying_matrices, ignore_existing_files=ignore_existing_files,
                                        progress_bar=progress_bar)
        matrix_list = self._get_matrix_list(underlying_matrices)
        self._store_matrix_list(matrix_list, shape)

    def _get_underlying_matrices(self):
        op_matrix_list = []
        for (domain, target) in itertools.product(self.domain.get_vs_list(), self.target.get_vs_list()):
            if self.operator_cls1.is_match(domain, target):
                op_matrix_list.append(self.operator_cls1(domain, target))
            if self.operator_cls2.is_match(domain, target):
                op_matrix_list.append(self.operator_cls2(domain, target))
        return op_matrix_list

    def _build_underlying_matrices(self, op_matrix_list, **kwargs):
        for op in op_matrix_list:
            op.build_matrix(**kwargs)

    def _get_matrix_list(self, underlying_matrices):
        matrixList = []
        for op in underlying_matrices:
            if not op.is_valid():
                continue
            domain_start_idx = self.domain.get_start_idx(op.get_domain())
            target_start_idx = self.target.get_start_idx(op.get_target())
            subMatrixList = op.get_matrix_list()
            for (i, j, v) in subMatrixList:
                matrixList.append((i + domain_start_idx, j + target_start_idx, v))
        matrixList.sort()
        return matrixList


class OperatorMatrixCollection(object):
    def __init__(self, sum_vector_space, op_matrix_list):
        self.sum_vector_space = sum_vector_space
        self.op_matrix_list = op_matrix_list
        self.info_tracker = None

    @abstractmethod
    def get_type(self):
        pass

    def __str__(self):
        return '<%s operator matrix collection on %s>' % (self.get_type(), str(self.sum_vector_space))

    def get_op_list(self):
        return self.op_matrix_list

    def get_vector_space(self):
        return self.sum_vector_space

    def sort(self, key='work_estimate'):
        if key == 'work_estimate':
            self.op_matrix_list.sort(key=operator.methodcaller('get_work_estimate'))
        elif key == 'size':
            self.op_matrix_list.sort(key=operator.methodcaller('get_sort_size'))
        elif key == 'entries':
            self.op_matrix_list.sort(key=operator.methodcaller('get_sort_entries'))
        else:
            raise ValueError("Invalid sort key. Options: 'work_estimate', 'size', 'entries'")

    def build_matrix(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        print(' ')
        print('Build matrices of %s' % str(self))
        if n_jobs > 1:
            info_tracker = False
            progress_bar = False
        if info_tracker:
            self.start_tracker()
        self.sort()
        Parallel.parallel(self._build_single_matrix, self.op_matrix_list, n_jobs=n_jobs,
                          ignore_existing_files=ignore_existing_files, info_tracker=info_tracker,
                          progress_bar=progress_bar)
        if info_tracker:
            self.stop_tracker()

    def _build_single_matrix(self, op, info_tracker=False, **kwargs):
        op.build_matrix(info_tracker=info_tracker, **kwargs)
        if info_tracker:
            self.update_tracker(op)

    def compute_rank(self, exact=False, n_primes=1, estimate=False, sort_key='size', ignore_existing_files=False,
                     n_jobs=1, info_tracker=False):
        print(' ')
        print('Compute ranks of %s' % str(self))
        if n_jobs > 1:
            info_tracker = False
        if info_tracker:
            self.start_tracker()
        self.sort(key=sort_key)
        Parallel.parallel(self._compute_single_rank, self.op_matrix_list, n_jobs=n_jobs, exact=exact, n_primes=n_primes,
                    estimate=estimate, ignore_existing_files=ignore_existing_files, info_tracker=info_tracker)
        if info_tracker:
            self.stop_tracker()

    def _compute_single_rank(self, op, info_tracker=False, **kwargs):
        op.compute_rank(**kwargs)
        if info_tracker:
            self.update_tracker(op)

    def set_tracker_parameters(self):
        try:
            param_names = self.get_vector_space().get_vs_list()[0].get_ordered_param_dict().keys()
        except IndexError:
            param_names = []
        parameter_list = param_names + OperatorMatrixProperties.names()
        self.info_tracker.set_parameter_list(parameter_list)

    def start_tracker(self):
        self.info_tracker = DisplayInfo.InfoTracker(str(self))
        self.set_tracker_parameters()
        op_info_dict = collections.OrderedDict()
        for op in self.op_matrix_list:
            op_info_dict.update({tuple(op.domain.get_ordered_param_dict().values()): op.get_properties().list()})
        self.info_tracker.update_data(op_info_dict)
        self.info_tracker.start()

    def update_tracker(self, op):
        op.update_properties()
        message = {tuple(op.domain.get_ordered_param_dict().values()): op.get_properties().list()}
        self.info_tracker.get_queue().put(message)

    def stop_tracker(self):
        self.info_tracker.stop()


class Differential(OperatorMatrixCollection):
    __metaclass__ = ABCMeta

    def __init__(self, sum_vector_space, op_matrix_list):
        super(Differential, self).__init__(sum_vector_space, op_matrix_list)

    @abstractmethod
    def get_cohomology_plot_path(self):
        pass

    @abstractmethod
    def get_cohomology_plot_parameter_order(self):
        pass

    def get_ordered_cohomology_param_range_dict(self):
        return self.sum_vector_space.get_ordered_param_range_dict()

    def __str__(self):
        return '<%s differential on %s>' % (self.get_type(), str(self.sum_vector_space))

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
            logger.info("Cannot compute cohomology: First build basis for %s " % str(opD.get_domain()))
            return None
        if dimV == 0:
            return 0
        if opD.is_valid():
            try:
                rankD = opD.get_matrix_rank()
            except SL.FileNotFoundError:
                logger.info("Cannot compute cohomology: Matrix rank not calculated for %s " % str(opD))
                return None
        else:
            rankD = 0
        if opDD.is_valid():
            try:
                rankDD = opDD.get_matrix_rank()
            except SL.FileNotFoundError:
                logger.info("Cannot compute cohomology: Matrix rank not calculated for %s " % str(opDD))
                return None
        else:
            rankDD = 0
        cohomologyDim = dimV - rankD - rankDD
        if cohomologyDim < 0:
            print("Negative cohomology dimension for %s" % str(opD.domain))
        return cohomologyDim

    # Computes the cohomology, i.e., ker(D)/im(DD)
    def get_general_cohomology_dim_dict(self):
        cohomology_dim = dict()
        for (opD, opDD) in itertools.permutations(self.op_matrix_list, 2):
            if Differential.is_match(opD, opDD):
                dim = Differential.cohomology_dim(opD, opDD)
                cohomology_dim.update({opD.domain: dim})
        return cohomology_dim

    def get_cohomology_dim(self):
        cohomology_dim = self.get_general_cohomology_dim_dict()
        dim_dict = dict()
        for vs in self.sum_vector_space.get_vs_list():
            dim_dict.update({vs.get_ordered_param_dict().get_value_tuple(): cohomology_dim.get(vs)})
        return dim_dict

    def square_zero_test(self, eps=Parameters.square_zero_test_eps):
        print(' ')
        print("Square zero test for %s:" % str(self))
        succ_l = 0  # number of pairs for which test was successful
        fail = []  # failed pairs
        triv_l = 0  # number of pairs for which test trivially succeeded because at least one operator is the empty matrix
        inc_l = 0  # number of pairs for which operator matrices are missing
        for (op1, op2) in itertools.permutations(self.op_matrix_list, 2):
            if Differential.is_match(op2, op1):
                # A composable pair is found

                pair = (op1, op2)
                res = self._square_zero_test_for_pair(pair, eps=eps)

                if res == 'triv':
                    triv_l += 1
                elif res == 'succ':
                    succ_l += 1
                elif res == 'inc':
                    inc_l += 1
                elif res == 'fail':
                    fail.append(pair)
                else:
                    raise ValueError('Undefined commutativity test result')

        fail_l = len(fail)
        print("trivial success: %d, success: %d, inconclusive: %d, failed: %d pairs" % (triv_l, succ_l, inc_l, fail_l))
        logger.warn("Square zero test for %s:" % str(self))
        logger.warn("trivial success: %d, success: %d, inconclusive: %d, failed: %d pairs" %
                    (triv_l, succ_l, inc_l, fail_l))
        for (op1, op2) in fail:
            logger.error("Square zero test for %s: failed for the pair %s, %s" % (str(self), str(op1), str(op2)))
        return (triv_l, succ_l, inc_l, fail_l)

    def _square_zero_test_for_pair(self, pair, eps=Parameters.square_zero_test_eps):
        (op1, op2) = pair
        if not (op1.is_valid() and op2.is_valid()):
            return 'triv'
        try:
            if op1.is_trivial() or op2.is_trivial():
                return 'triv'
            M1 = op1.get_matrix()
            M2 = op2.get_matrix()
        except SL.FileNotFoundError:
            logger.info("Cannot test square zero: "
                        "Operator matrix not built for %s or %s" % (str(op1), str(op2)))
            return 'inc'
        if SH.matrix_norm(M2 * M1) < eps:
            return 'succ'
        else:
            return 'fail'

    def plot_cohomology_dim(self):
        dim_dict = self.get_cohomology_dim()
        plot_path = self.get_cohomology_plot_path()
        param_order = self.get_cohomology_plot_parameter_order()
        ordered_param_range_dict = self.get_ordered_cohomology_param_range_dict()
        PlotCohomology.plot_array(dim_dict, ordered_param_range_dict, plot_path, param_order)


class BiDifferential(Differential):
    __metaclass__ = ABCMeta

    def __init__(self, graded_sum_vs, operator_cls1, operator_cls2, bi_op_matrix_cls):
        self.graded_sum_vs = graded_sum_vs
        self.operator_cls1 = operator_cls1
        self.operator_cls2 = operator_cls2
        self.bi_op_matrix_cls = bi_op_matrix_cls
        op_matrix_list = self.bi_op_matrix_cls.generate_op_matrix_list(graded_sum_vs, self.operator_cls1, self.operator_cls2)
        super(BiDifferential, self).__init__(self.graded_sum_vs, op_matrix_list)
