"""Provide abstract classes for operator matrices, graph operators, operator matrix collections, differentials,
and bi operator matrices to be used in bicomplexes"""

__all__ = ['OperatorMatrixProperties', 'OperatorMatrix', 'Operator', 'GraphOperator', 'BiOperatorMatrix',
           'OperatorMatrixCollection', 'Differential']

from abc import ABCMeta, abstractmethod
import itertools
from tqdm import tqdm
import collections
import scipy.sparse as sparse
from sage.all import *
import Log
import StoreLoad
import Parallel
import Shared
import Parameters
import PlotCohomology
import DisplayInfo
import RheinfallInterface
import LinboxInterface
import GraphVectorSpace

logger = Log.logger.getChild('graph_operator')


class OperatorMatrixProperties(object):
    """Properties of an operator matrix.

    Attributes:
        - valid (bool): Validity of the operator matrix.
        - shape (tuple(int, int)): Shape of the operator matrix.
        - entries (int): Number of nonzero entries of the sparse matrix.
        - rank (int): Rank of the operator matrix.
    """

    def __init__(self):
        """Initialize the matrix properties with None."""
        self.valid = None
        self.shape = None
        self.entries = None
        self.rank = None

    @classmethod
    def names(cls):
        """Return a list of the names of the matrix properties.

        :return: Names of the matrix properties.
        :rtype: list(str)
        """
        return ['valid', 'shape', 'entries', 'rank']

    @staticmethod
    def sort_variables():
        """Return a list of the matrix properties used as sort keys for operator matrices.

        :return: Names of the matrix properties, used as sort keys for operator matrices.
        :rtype: list(str)
        """
        return ['valid', 'entries']

    def list(self):
        """Return a list of the matrix properties.

        :return: Matrix properties.
        :rtype: list(int)
        """
        return [self.valid, self.shape, self.entries, self.rank]


class OperatorMatrix(object):
    """Operator matrix class.

    Abstract class defining the interface for an operator matrix.

    Attributes:
        - domain (GraphVectorSpace.VectorSpace): Domain vector space.
        - target (GraphVectorSpace.VectorSpace): Target vector space.
        - properties (OperatorMatrixProperties): Operator matrix properties object, containing information like
            validity, shape, number of nonzero entries, and rank.
    """
    __metaclass__ = ABCMeta

    # Data type to store the matrix using the SMS (http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html) format.
    data_type = "M"

    def __init__(self, domain, target):
        """Initialize domain and target graph vector spaces and initialize the operator matrix properties with None.

        :param domain: Domain vector space.
        :type domain: GraphVectorSpace.VectorSpace
        :param target: Target vector space.
        :type target: GraphVectorSpace.VectorSpace
        :raise: ValueError: Raised if domain and target don't match to build the operator matrix.
        """
        if not self.is_match(domain, target):
            raise ValueError("Domain %s and target %s don't match to build the operator matrix %s"
                             % (str(domain), str(target), str(self)))
        self.domain = domain
        self.target = target
        self.properties = OperatorMatrixProperties()
        # set is_pseudo_matrix=True if the matrix actually stored is not operating between the vector spaces.
        self.is_pseudo_matrix = False

    def get_domain(self):
        """Return the domain graph vector space of the operator matrix.

        :return: Domain vector space.
        :rtype: GraphVectorSpace.VectorSpace
        """
        return self.domain

    def get_target(self):
        """Return the target graph vector space of the operator matrix.

        :return: Target vector space.
        :rtype: GraphVectorSpace.VectorSpace
        """
        return self.target

    @abstractmethod
    def __str__(self):
        """Return a unique description of the graph operator matrix.

        :return: Unique description of the graph operator matrix.
        :rtype: str
        """
        pass

    @abstractmethod
    def get_matrix_file_path(self):
        """Return the path for the operator matrix file.

        :return: Path to the operator matrix file.
        :rtype: path
        """
        pass

    @abstractmethod
    def get_rank_file_path(self):
        """Return the path for the matrix rank file.

        :return: Path to the matrix rank file.
        :rtype: path
        """
        pass

    def get_ref_matrix_file_path(self):
        """Return the path for the reference operator matrix file.

        Refer to reference data (if available) for testing.

        :return: Path to the reference operator matrix file.
        :rtype: path
        """
        pass

    def get_ref_rank_file_path(self):
        """Return the path for the reference matrix rank file.

        Refer to reference data (if available) for testing.

        :return: Path to the reference matrix rank file.
        :rtype: path
        """
        pass

    def is_valid(self):
        """Return the validity of the parameter combination for the operator matrix.

        :return: True if the domain and target vector spaces are valid, False otherwise.
        :rtype: bool
        """
        return self.domain.is_valid() and self.target.is_valid()

    @abstractmethod
    def get_work_estimate(self):
        """Estimate the work needed to build the operator matrix.

        Arbitrary units. Used to schedule the order of building the operator matrices.

        :return: Estimate the work to build the operator matrix. Arbitrary units.
        :rtype: int
        """
        pass

    @abstractmethod
    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        """Build the operator matrix and write it to the matrix file if the graph operator is valid.

        The operator matrix is stored in the SMS format to the matrix file. The header line contains the shape of the
        matrix (nrows = domain dimension, ncols = target dimension) and the data type of the SMS format. In the file the
        matrix entries are listed as (domain index, target index, value).

        :param ignore_existing_files: Option to ignore an existing matrix file. Ignore an existing file and
            rebuild the operator matrix if True, otherwise skip rebuilding the matrix file if there exists a
            matrix file already (Default: False).
        :type ignore_existing_files: bool
        :param skip_if_no_basis: Skip building the matrix if the domain or target basis is not built
            (Default: True). If False and a basis file is missing an error is raised.
        :type skip_if_no_basis: bool
        :param progress_bar: Option to show a progress bar (Default: False).
        :type progress_bar: bool
        :param kwargs: Accept further keyword arguments without influence.
        :raise StoreLoad.FileNotFoundError: Raised if skip_if_no_basis is True and the domain or target basis file is
            not found.

        .. seealso:: http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
        """
        pass

    @staticmethod
    @abstractmethod
    def is_match(domain, target):
        """Return whether domain and target vector space match to build an operator matrix.

        :param domain: Potential domain vector space for an operator matrix.
        :type domain: GraphVectorSpace.VectorSpace
        :param target: Potential target vector space for an operator matrix.
        :type targetGraphVectorSpace.VectorSpace
        :return:True if domain and target match to build an operator matrix.
        :rtype: bool
        """
        pass

    def exists_matrix_file(self):
        """Return whether there exists a matrix file.

        :return: True if a matrix file is found.
        :rtype: bool
        """
        return os.path.isfile(self.get_matrix_file_path())

    def exists_rank_file(self):
        """Return whether there exists a rank file.

        :return: True if a rank file is found.
        :rtype: bool
        """
        return os.path.isfile(self.get_rank_file_path())

    def delete_matrix_file(self):
        """Delete the matrix file."""
        if os.path.isfile(self.get_matrix_file_path()):
            os.remove(self.get_matrix_file_path())

    def delete_rank_file(self):
        """Delete the rank file."""
        if os.path.isfile(self.get_rank_file_path()):
            os.remove(self.get_rank_file_path())

    def _store_matrix_list(self, matrix_list, shape, data_type=data_type):
        """Store the operator matrix in SMS format to the matrix file.

        The header line contains the shape of the matrix (nrows = domain dimension, ncols = target dimension) and
        the data type of the SMS format. In the file the matrix entries are listed as (domain index, target index, value).

        :param matrix_list: List of matrix entries in the form (domain index, target index, value).
            The list entries must be ordered lexicographically.
        :type matrix_list: list(tuple(int, int, int))
        :param shape: Tuple containing the matrix shape = (nrows = domain dimension, ncols = target dimension).
        :type shape: tuple(int, int)
        :param data_type: data_type for the SMS format.
        :type data_type: str

        ..seealso:: http://ljk.imag.fr/membres/Jean-Guillaume.Dumas/simc.html
        """
        (d, t) = shape
        stringList = []
        stringList.append("%d %d %s" % (d, t, data_type))
        for (i, j, v) in matrix_list:
            stringList.append("%d %d %d" % (i + 1, j + 1, v))
        stringList.append("0 0 0")
        StoreLoad.store_string_list(stringList, self.get_matrix_file_path())

    def _load_matrix_list(self):
        """Load the operator matrix in the SMS format from the matrix file.

        Return the matrix list, i.e. a list of the non-zero matrix entries, and the matrix shape, i.e. a tuple with
        the number of rows and columns.
        :return: : (matrix_list = list((domain index, target index, value), shape = (domain dimension, target dimension))
        :rtype: tuple(list(tuple(int, int, int)), tuple(int, int))
        :raise StoreLoad.FileNotFoundError: Raised if the matrix file cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        if not self.exists_matrix_file():
            raise StoreLoad.FileNotFoundError(
                "Cannot load matrix, No matrix file found for %s: " % str(self))
        stringList = StoreLoad.load_string_list(self.get_matrix_file_path())
        (d, t, data_type) = stringList.pop(0).split(" ")
        shape = (d, t) = (int(d), int(t))
        if (not self.is_pseudo_matrix) and (d != self.domain.get_dimension() or t != self.target.get_dimension()):
            raise ValueError("%s: Shape of matrix doesn't correspond to the vector space dimensions"
                             % str(self.get_matrix_file_path()))
        tail = map(int, stringList.pop().split(" "))
        if not list(tail) == [0, 0, 0]:
            raise ValueError("%s: End line missing or matrix not correctly read from file"
                             % str(self.get_matrix_file_path()))
        matrix_list = []
        for line in stringList:
            (i, j, v) = map(int, line.split(" "))
            if i < 1 or j < 1:
                raise ValueError("%s: Invalid matrix index: %d %d" %
                                 (str(self.get_matrix_file_path()), i, j))
            if i > d or j > t:
                raise ValueError("%s: Invalid matrix index outside matrix size:"
                                 " %d %d" % (str(self.get_matrix_file_path()), i, j))
            matrix_list.append((i - 1, j - 1, v))
        return (matrix_list, shape)

    def get_matrix_list(self):
        """Return the list of nonzero matrix entries.

        :return: List of matrix entries in the form (domain index, target index, value).
        :rtype: list(tuple(int, int, int))
        :raise StoreLoad.FileNotFoundError: Raised if the matrix file cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        (matrix_list, shape) = self._load_matrix_list()
        return matrix_list

    def get_shifted_matrix_list(self, domain_start, target_start):
        """Return the list of nonzero matrix entries with indices shifted by the domain and target start indices.

        The shifted matrix list is used to build the operator matrix of bioperators in bicomplexes.

        :param domain_start: Domain start index.
        :type domain_start: int
        :param target_start: Target start index.
        :type target_start: int
        :return: List of matrix entries in the form
            (domain index + domain start index, target index + target start index, value).
        :rtype: : list(tuple(int, int, int))
        :raise StoreLoad.FileNotFoundError: Raised if the matrix file cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        matrix_list = self.get_matrix_list()
        shifted_matrix_list = []
        for (i, j, v) in matrix_list:
            shifted_matrix_list.append((i + domain_start, j + target_start, v))
        return shifted_matrix_list

    def get_matrix_shape(self):
        """Return the matrix shape.

        :return: Matrix shape = (nrows = target dimension, ncolumns = domain dimension).
        :rtype: tuple(int, int)
        :raise StoreLoad.FileNotFoundError: Raised if there is neither a matrix file nor the basis files of the domain and
            target vector spaces.
        """
        try:
            header = StoreLoad.load_line(self.get_matrix_file_path())
            (d, t, data_type) = header.split(" ")
            (d, t) = (int(d), int(t))
        except StoreLoad.FileNotFoundError:
            try:
                d = self.domain.get_dimension()
                t = self.target.get_dimension()
            except StoreLoad.FileNotFoundError:
                raise StoreLoad.FileNotFoundError("Matrix shape of %s unknown: "
                                                  "Build matrix or domain and target basis first" % str(self))
        return (t, d)

    def get_matrix_shape_entries(self):
        """Return the matrix shape and the number of non-zero matrix entries.

        :return: (matrix shape, entries).
                Matrix shape = (nrows = target dimension, ncolumns = domain dimension) and number of non-zero matrix
                entries.
        :rtype: tuple(tuple(int, int), int)
        :raise StoreLoad.FileNotFoundError: Raised if the matrix file cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        try:
            (matrixList, shape) = self._load_matrix_list()
            (d, t) = shape
            return ((t, d), len(matrixList))
        except StoreLoad.FileNotFoundError:
            raise StoreLoad.FileNotFoundError(
                "Matrix shape and entries unknown for %s: No matrix file" % str(self))

    def get_matrix_entries(self):
        """Return the number of non-zero matrix entries.

        :return: int: Number of non-zero matrix entries.
        :rtype: int
        :raise StoreLoad.FileNotFoundError: Raised if the matrix is valid and the matrix file cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        if not self.is_valid():
            return 0
        (shape, entries) = self.get_matrix_shape_entries()
        return entries

    def is_trivial(self):
        """Return whether the matrix is trivial, i.e. is not valid, has zero dimension or hasn't any non-zero entries

        :return: True if the matrix is trivial (not valid, zero dimension or no non-zero entries).
        :rtype: bool
        :raise StoreLoad.FileNotFoundError: Raised if the matrix is valid and the matrix file or the basis files of domain and
                target cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        if not self.is_valid():
            return True
        (t, d) = self.get_matrix_shape()
        if t == 0 or d == 0:
            return True
        if self.get_matrix_entries() == 0:
            return True
        return False

    def get_matrix_transposed(self):
        """Return the transposed operator matrix as sparse sage matrix over Z.

        :return: Transposed operator matrix with shape (domain dimension, target dimension).
        :rtype: Matrix_sparse
        :raise StoreLoad.FileNotFoundError: Raised if the matrix is valid and the matrix file cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        if not self.is_valid():
            logger.warning("Zero matrix: %s is not valid" % str(self))
            (d, t) = (self.domain.get_dimension(), self.target.get_dimension())
            entries_list = []
        else:
            (entries_list, shape) = self._load_matrix_list()
            (d, t) = shape
        M = matrix(ZZ, d, t, sparse=True)
        for (i, j, v) in entries_list:
            M.add_to_entry(i, j, v)
        return M

    def get_matrix(self):
        """Return the operator matrix as sparse sage matrix over Z.

        :return: Operator matrix with shape (target dimension, domain dimension).
        :rtype: Matrix_sparse
        :raise StoreLoad.FileNotFoundError: Raised if the matrix is valid and the matrix file cannot be found or if
            the matrix is not valid and the domain and target basis files cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        M = self.get_matrix_transposed().transpose()
        if (not self.is_pseudo_matrix) and (M.ncols() != self.get_domain().get_dimension() or M.nrows() != self.get_target().get_dimension()):
            raise ValueError(
                "Matrix shape doesn't match the dimension of the domain or the target for " + str(self))
        return M

    def get_matrix_scipy_transposed(self):
        """Return the transposed operator matrix as sparse csc matrix.

        :return: Transposed operator matrix with shape (domain dimension, target dimension).
        :rtype: scipy.sparse.csc_martix
        :raise StoreLoad.FileNotFoundError: Raised if the matrix is valid and the matrix file cannot be found or if
            the matrix is not valid and the domain and target basis files cannot be found.
        :raise ValueError: Raised in the following cases:
                The shape of the matrix doesn't correspond to the dimensions of the domain or target vector space.
                End line is missing.
                Non-positive matrix indices.
                Matrix indices outside matrix shape.
        """
        data = []
        row_ind = []
        col_ind = []
        if not self.is_valid():
            logger.warning("Zero matrix: %s is not valid" % str(self))
            shape = (self.domain.get_dimension(), self.target.get_dimension())
        else:
            (entries_list, shape) = self._load_matrix_list()
            for (r, c, d) in entries_list:
                row_ind.append(r)
                col_ind.append(c)
                data.append(d)
        M = sparse.csc_matrix((data, (row_ind, col_ind)),
                              shape=shape, dtype='d')
        return M

    def compute_rank(self, sage=None, linbox=None, rheinfall=None, ignore_existing_files=False, skip_if_no_matrix=True):
        """Compute the rank of the operator matrix.

        Compute the rank of the operator matrix and stores it in the rank file. The rank can be determined with
        different modes:
            - Use sage to determine the rank over the integers or over a finite field, i.e. with all
              calculations modulo a prime number.
            - Use linbox to determine the rank over the rational numbers or over a finite field, i.e. with all
              calculations modulo a prime number.
            - Use rheinfall to determine the rank over the integers, the rational numbers or modulo 64 bit integers.
        The prime number is set in the module Parameters.

        :param sage: Use sage to compute the rank. Options: 'integer' (exact rank over the integers),
            'mod' (rank over a finite field, i.e. calculations modulo a prime number (Default: None).
        :type sage: str or list(str)
        :param linbox: Use linbox to compute the rank. Options: 'rational' (exact rank over the rational numbers),
            'mod' (rank over a finite field, i.e. calculations modulo a prime number (Default: None).
        :type linbox: str or list(str)
        :param rheinfall: Use rhainfall to compute the rank. Options: 'int64', 'mpq', 'mpz' (Default: None).
        :type rheinfall: str or list(str)
        :param ignore_existing_files: Option to ignore an existing rank file. Ignore existing file and
            recompute the rank if True, otherwise skip recomputing the rank if there exists already a
            rank file (Default: False).
        :type ignore_existing_files: bool
        :raise StoreLoad.FileNotFoundError: Raised if the matrix file cannot be found and skip_if_no_matrix = False.

        .. seealso:: - http://www.linalg.org/
                    - https://github.com/linbox-team/linbox/blob/master/examples/rank.C
                    - https://github.com/riccardomurri/rheinfall/blob/master/src.c%2B%2B/examples/rank.cpp
        """
        if sage is None and linbox is None and rheinfall is None:
            raise ValueError("compute_rank: At least one rank computation method needs to be specified.")

        if not self.is_valid():
            return
        if not ignore_existing_files and self.exists_rank_file():
            return
        elif self.exists_rank_file():
            self.delete_rank_file()
        print('Compute matrix rank: Domain: ' +
              str(self.domain.get_ordered_param_dict()))
        try:
            rank_dict = self._compute_rank(
                sage=sage, linbox=linbox, rheinfall=rheinfall)
        except StoreLoad.FileNotFoundError as error:
            if skip_if_no_matrix:
                logger.info(
                    "Skip computing rank of %s, since matrix is not built" % str(self))
                # print("Exception: "+str(error))
                return
            else:
                raise error
        self._store_rank_dict(rank_dict)

    def _compute_rank(self, sage=None, linbox=None, rheinfall=None, prime=Parameters.prime):
        if type(sage) == str:
            sage = [sage]
        if type(linbox) == str:
            linbox = [linbox]
        if type(rheinfall) == str:
            rheinfall = [rheinfall]
        if self.is_trivial() or self.get_matrix_entries() == 0:
            rank_dict = {'exact': 0}
        else:
            rank_dict = {}
            try:
                if sage is not None:
                    if not set(sage) <= Parameters.sage_rank_options:
                        raise ValueError(
                            "Options for rank computations with sage: " + str(Parameters.sage_rank_options))
                    for option in sage:
                        if option == 'integer':
                            M = self.get_matrix_transposed()
                            rank_exact = M.rank()
                            rank_dict.update({'sage_integer': rank_exact})
                        if option == 'mod':
                            M = self.get_matrix_transposed()
                            M.change_ring(GF(prime))
                            rank_mod_p = M.rank()
                            info = 'sage_mod_%d' % prime
                            rank_dict.update({info: rank_mod_p})
                if linbox is not None:
                    if not set(linbox) <= LinboxInterface.linbox_options:
                        raise ValueError("Options for rank computations with linbox: " +
                                         str(LinboxInterface.linbox_options))
                    for option in linbox:
                        rank_linbox = LinboxInterface.rank(
                            option, self.get_matrix_file_path(), prime=prime)
                        info = "linbox_%s" % option if option == "rational" else "linbox_%s_%d" % (
                            option, prime)
                        rank_dict.update({info: rank_linbox})
                if rheinfall is not None:
                    if not set(rheinfall) <= RheinfallInterface.rheinfall_options:
                        raise ValueError("Options for rank computations with rheinfall: " +
                                         str(RheinfallInterface.rheinfall_options))
                    for option in rheinfall:
                        rank_rheinfall = RheinfallInterface.rank(
                            option, self.get_matrix_file_path())
                        if rank_rheinfall == 0:
                            return self._compute_rank(sage='integer')
                        info = "rheinfall_" + option
                        rank_dict.update({info: rank_rheinfall})
            except StoreLoad.FileNotFoundError:
                raise StoreLoad.FileNotFoundError(
                    "Cannot compute rank of %s: First build operator matrix" % str(self))
        return rank_dict

    def exists_exact_rank(self):
        """Determines whether there has been a rank computed that is exact, i.e., 
        over rationals."""
        if not self.is_valid():
            return True
        if not self.exists_rank_file():
            return False

        rank_dict = self._load_rank_dict()
        return "sage_integer" in rank_dict or "exact" in rank_dict or "linbox_rational" in rank_dict

    def _store_rank_dict(self, update_rank_dict):
        try:
            rank_dict = self._load_rank_dict()
        except StoreLoad.FileNotFoundError:
            rank_dict = dict()
        rank_dict.update(update_rank_dict)
        rank_list = [str(rank) + ' ' + mode for (mode, rank)
                     in rank_dict.items()]
        StoreLoad.store_string_list(rank_list, self.get_rank_file_path())

    def _load_rank_dict(self):
        if not self.is_valid():
            return {'exact': 0}
        try:
            rank_list = StoreLoad.load_string_list(self.get_rank_file_path())
        except StoreLoad.FileNotFoundError:
            raise StoreLoad.FileNotFoundError(
                "Cannot load matrix rank, No rank file found for %s: " % str(self))
        rank_dict = dict()
        for line in rank_list:
            (rank, mode) = line.split(" ")
            rank_dict.update({mode: int(rank)})
        return rank_dict

    def get_matrix_rank(self):
        """Return the matrix rank.

        :return: Matrix rank.
        :rtype: int
        :raise StoreLoad.FileNotFoundError: Raised if the rank file is not found.
        """
        if (not self.is_valid()) or self.domain.get_dimension() == 0 or self.target.get_dimension() == 0:
            return 0

        rank_dict = self._load_rank_dict()
        ranks = list(rank_dict.values())
        if len(ranks) == 0:
            raise ValueError(
                "No matrix rank stored in rank file for " + str(self))
        if len(set(ranks)) != 1:
            raise ValueError(
                "Matrix ranks computed with different methods are not equal for " + str(self))
        return ranks[0]

    def get_sort_size(self):
        """Return the min(nrows, ncolumns) to be used as a sort key. If the matrix shape is unknown the constant
        Parameters.max_sort_value is returned.

        :return: min(nrows, ncolumns) or a constant to be used as a sort key
        :rtype: int
        """
        try:
            sort_size = min(self.get_matrix_shape())
        except StoreLoad.FileNotFoundError:
            sort_size = Parameters.max_sort_value
        return sort_size

    def get_sort_entries(self):
        """Return the number of nonzero entries to be used as a sort key. If the number of entries is unknown the
        constant Parameters.max_sort_value is returned.

        :return: Number of non-zero entries or a constant to be used as a sort key
        :rtype: int
        """
        try:
            sort_entries = self.get_matrix_entries()
        except StoreLoad.FileNotFoundError:
            sort_entries = Parameters.max_sort_value
        return sort_entries

    def update_properties(self):
        """Update the vector space properties by reading from the matrix and rank files if available.
        """
        self.properties.valid = self.is_valid()
        try:
            self.properties.shape = self.get_matrix_shape()
        except StoreLoad.FileNotFoundError:
            pass
        try:
            self.properties.entries = self.get_matrix_entries()
        except StoreLoad.FileNotFoundError:
            pass
        try:
            self.properties.rank = self.get_matrix_rank()
        except StoreLoad.FileNotFoundError:
            pass

    def get_properties(self):
        """Return the operator matrix properties.

        :return: Operator matrix properties.
        :rtype: OperatorMatrixProperties
        """
        return self.properties


class Operator(object):
    """Interface for operators."""
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_type(self):
        """Return a unique description of the operator type.

        :return: Unique description of the operator type.
        :rtype: str
            """
        pass

    @abstractmethod
    def operate_on(self, G):
        """For G a sage graph returns a list of pairs (GG, factor), such that (operator)(G) = sum(factor * GG).

        :param G: Graph on which the operator is applied.
        :type G: Graph
        :return: List of tuples (GG, factor), such that (operator)(G) = sum(factor * GG)
        :rtype: list(tuple(Graph, factor))
        """
        pass


class GraphOperator(Operator, OperatorMatrix):
    """Graph oerator.

    Inherits from operator matrix and additionally implements the operator interface, i.e. it acts on graphs as
    described in the method operate on. Build the operator matrix by calling the method build_matrix.
    """
    __metaclass__ = ABCMeta

    def __init__(self, domain, target):
        super(GraphOperator, self).__init__(domain, target)

    @classmethod
    def generate_op_matrix_list(cls, sum_vector_space):
        """Return a list of all possible graph operators of this type with domain and target being sub vector spaces
        of the sum vector space.

        :param sum_vector_space: Sum vector space of graphs on which the graph
            operator is based.
        :type sum_vector_space: GraphVectorSpace.SumVectorSpace
        :return: List of all possible graph operators with domain and target in the sum vector space.
        :rtype: list(GraphOperator)
        """
        op_matrix_list = []
        for (domain, target) in itertools.permutations(sum_vector_space.get_vs_list(), 2):
            if cls.is_match(domain, target):
                op_matrix_list.append(cls(domain, target))
        return op_matrix_list

    def __str__(self):
        """Return a unique description of the graph operator.

        :return: Unique description of the graph operator.
        :rtype: str
        """
        return '<%s graph operator, domain: %s>' % (self.get_type(), str(self.domain))

    def operate_on_list(self, graph_sgn_list):
        """Operate on the output of the operate on method of the graph operator.

        This method is intended to act with the graph operator twice on a graph.

        :param graph_sgn_list: Output of the operate_on method of the same graph operator.
        :type graph_sgn_list: list(tuple(Graph, factor))
        :return: List of tuples (GG, factor), such that (operator)(graph_sgn_list) = sum(factor * GG)
        :rtype: list(tuple(Graph, factor))
        """
        image_dict = dict()
        for (G1, sgn1) in graph_sgn_list:
            G1_image = self.operate_on(G1)
            for (G2, sgn2) in G1_image:
                (G2_g6, sgn_p) = self.target.graph_to_canon_g6(G2)
                v = image_dict.get(G2_g6)
                new_v = sgn1 * sgn2 * sgn_p
                if v is not None:
                    new_v += v
                if new_v:
                    image_dict.update({G2_g6: new_v})
                elif v is not None:
                    image_dict.pop(G2_g6)
        return image_dict.items()

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=False, **kwargs):
        if not self.is_valid():
            return
        if (not ignore_existing_files) and self.exists_matrix_file():
            return
        try:
            domain_basis = self.domain.get_basis()
        except StoreLoad.FileNotFoundError:
            if not skip_if_no_basis:
                raise StoreLoad.FileNotFoundError("Cannot build operator matrix of %s: "
                                                  "First build basis of the domain %s" % (str(self), str(self.domain)))
            else:
                logger.info("Skip building operator matrix of %s "
                            "since basis of the domain %s is not built" % (str(self), str(self.domain)))
                return
        try:
            target_basis6 = self.target.get_basis_g6()
        except StoreLoad.FileNotFoundError:
            if not skip_if_no_basis:
                raise StoreLoad.FileNotFoundError("Cannot build operator matrix of %s: "
                                                  "First build basis of the target %s" % (str(self), str(self.target)))
            else:
                logger.info("Skip building operator matrix of %s "
                            "since basis of the target %s is not built" % (str(self), str(self.target)))
                return

        shape = (d, t) = (self.domain.get_dimension(),
                          self.target.get_dimension())
        if d == 0 or t == 0:
            self._store_matrix_list([], shape)
            return

        lookup = {G6: j for (j, G6) in enumerate(target_basis6)}
        desc = 'Build matrix of %s operator: Domain: %s' % (
            str(self.get_type()), str(self.domain.get_ordered_param_dict()))

        # if not progress_bar:
        print(desc)
        # list_of_lists = []
        matrix_list = []
        for domain_basis_element in tqdm(enumerate(domain_basis), total=d, desc=desc, disable=(not progress_bar)):
            # for domain_basis_element in list(enumerate(domain_basis)):
            # list_of_lists.append(self._generate_matrix_list(
            # domain_basis_element, lookup))
            matrix_list.extend(self._generate_matrix_list(
                domain_basis_element, lookup))
        # matrix_list = list(itertools.chain.from_iterable(list_of_lists))
        matrix_list.sort()
        self._store_matrix_list(matrix_list, shape)

    def _generate_matrix_list(self, domain_basis_element, lookup):
        (domain_index, G) = domain_basis_element
        image_list = self.operate_on(G)
        canon_images = dict()
        for (GG, prefactor) in image_list:
            (GGcanon6, sgn1) = self.target.graph_to_canon_g6(GG)
            sgn0 = canon_images.get(GGcanon6)
            sgn0 = sgn0 if sgn0 is not None else 0
            canon_images.update({GGcanon6: (sgn0 + sgn1 * prefactor)})
        matrix_list = []
        for (image, factor) in canon_images.items():
            if factor:
                target_index = lookup.get(image)
                if target_index is not None:
                    matrix_list.append((domain_index, target_index, factor))
        return matrix_list


class BiOperatorMatrix(OperatorMatrix):
    """Bi operator matrix to be used as operator matrix in bicomplexes."""
    __metaclass__ = ABCMeta

    @abstractmethod
    def __init__(self, domain, target, operator_cls1, operator_cls2):
        """Initialize domain and target degree slice and the two operator classes composing the bi operator.

        :param domain: Domain degree slice.
        :type domain: GraphVectorSpace.DegSlice
        :param target: Target degree slice.
        :type target: GraphVectorSpace.DegSlice
        :param operator_cls1: First operator class to compose the bi operator.
        :type operator_cls1: OperatorMatrix
        :param operator_cls2: Second operator class to compose the bi operator.
        :type operator_cls2: OperatorMatrix
        """
        super(BiOperatorMatrix, self).__init__(domain, target)
        self.operator_cls1 = operator_cls1
        self.operator_cls2 = operator_cls2

    @classmethod
    def generate_op_matrix_list(cls, graded_sum_vs):
        """Return a list of all possible bi operator matrices of this type with domain and target being degree slices
        of the graded sum vector space.

        :param graded_sum_vs: Graded sum vector space composed of degree slices.
        :type graded_sum_vs: GraphVectorSpace.SumVectorSpace
        :return: List of all possible bi operator matrices with domain and target being degree slices of the
            graded sum vector space.
        :rtype: list(BiOperatorMatrix)
        """
        graded_sum_vs_list = graded_sum_vs.get_vs_list()
        bi_op_matrix_list = []
        for (domain, target) in itertools.permutations(graded_sum_vs_list, 2):
            if cls.is_match(domain, target):
                bi_op_matrix_list.append(cls(domain, target))
        return bi_op_matrix_list

    def __str__(self):
        """Return a unique description of the bi operator matrix.

        :return: Unique description of the bi operator matrix.
        :rtype: str
        """
        return '<Bi operator matrix on domain: %s>' % str(self.domain)

    def is_valid(self):
        # return True
        return self.domain.is_valid() and self.target.is_valid()

    def get_work_estimate(self):
        """Estimate the work needed to build the bi operator matrix by the product of the dimensions of domain and target.

        Used to schedule the order of building the operator matrices.

        :return int: Estimate the work to build the operator matrix.
        :rtype: int
        """

        return self.domain.get_dimension() * self.target.get_dimension()

    def build_matrix(self, ignore_existing_files=False, progress_bar=False, **kwargs):
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
                matrixList.append(
                    (i + domain_start_idx, j + target_start_idx, v))
        matrixList.sort()
        return matrixList


class OperatorMatrixCollection(object):
    """Graph operator on the direct sum of graph vector spaces.

    Collection of operator matrices composing an operator on a sum vector space.

    Attributes:
        - sum_vector_space (GraphVectorSpace.SumVectorSpace): Underlying sum vector space.
        - op_matrix_list (list(OperatorMatrix)): List of operator matrices composing the operator.
        - info_tracker (DisplayInfo.InfoTracker): Tracker for information about the operator matrices in op_matrix_list.
            Tracker is only active if the different operator matrices are not built in parallel.
        """

    def __init__(self, sum_vector_space, op_matrix_list):
        """Initialize the underlying sum vector space and the list of operator matrices composing the operator.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: GraphVectorSpace.SumVectorSpace
        :param op_matrix_list: List of operator matrices composing the operator
        :type op_matrix_list: list(OperatorMatrix)
        """
        self.sum_vector_space = sum_vector_space
        self.op_matrix_list = op_matrix_list
        self.info_tracker = None

    @abstractmethod
    def get_type(self):
        """Return a unique description of the operator type.

        :return: Unique description of the operator type.
        :rtype: str

        :Example:
        'contract edges'
        """
        pass

    @abstractmethod
    def get_info_plot_path(self):
        pass

    def __str__(self):
        """Return a unique description of the operator.

        :return: Unique description of the operator.
        :rtype: str
        """
        return '<%s operator matrix collection on %s>' % (self.get_type(), str(self.sum_vector_space))

    def get_op_list(self):
        """Return the operator matrix list composing the operator.

        :return: List of operator matrices composing the operator.
        :rtype: list(OperatorMatrix)
        """
        return self.op_matrix_list

    def get_vector_space(self):
        """Return the underlying vector space.

        :return: Underlying direct sum of graph vector spaces.
        :rtype: GraphVectorSpace.SumVectorSpace
        """
        return self.sum_vector_space

    def sort(self, key='work_estimate'):
        """Sort the list of sub vector spaces.

        Possible sort keys: 'work_estimate', 'size', 'entries'.

        :param key: Sort key. Options: 'work_estimate', 'size', 'entries'.
        :type key: str
        :raise  ValueError: If an unknown sort key is given.
        """
        if key == 'work_estimate':
            self.op_matrix_list.sort(
                key=operator.methodcaller('get_work_estimate'))
        elif key == 'size':
            self.op_matrix_list.sort(
                key=operator.methodcaller('get_sort_size'))
        elif key == 'entries':
            self.op_matrix_list.sort(
                key=operator.methodcaller('get_sort_entries'))
        else:
            raise ValueError(
                "Invalid sort key. Options: 'work_estimate', 'size', 'entries'")

    def build_matrix(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        """Build the operator matrices of the operator.

        :param ignore_existing_files: Option to ignore  existing matrix files. Ignore existing files and
            rebuild the operator matrices if True, otherwise skip rebuilding a matrix file if there exists already a
            matrix file (Default: False).
        :type ignore_existing_files: bool
        :param n_jobs: Option to build different matrices in parallel using
                n_jobs parallel processes (Default: 1).
        :type n_jobs: int
        :param progress_bar: Option to show a progress bar (Default: False).
        :type progress_bar: bool
            Only active if different matrices are not built in parallel.
        :param info_tracker: Option to plot information about the sub vector spaces in a web page
            (Default: False). Only active if different matrices are not built in parallel.
        :type info_tracker: bool
        """
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
        info = info_tracker if Parameters.second_info else False
        op.build_matrix(info_tracker=info, **kwargs)
        if info_tracker:
            self.update_tracker(op)

    def compute_rank(self, sage=None, linbox=None, rheinfall=None, sort_key='size', ignore_existing_files=False,
                     n_jobs=1, info_tracker=False):
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

        print(' ')
        print('Compute ranks of %s' % str(self))
        if n_jobs > 1:
            info_tracker = False
        if info_tracker:
            self.start_tracker()
        self.sort(key=sort_key)
        Parallel.parallel(self._compute_single_rank, self.op_matrix_list, n_jobs=n_jobs, sage=sage, linbox=linbox,
                          rheinfall=rheinfall, ignore_existing_files=ignore_existing_files, info_tracker=info_tracker)
        if info_tracker:
            self.stop_tracker()

    def _compute_single_rank(self, op, info_tracker=False, **kwargs):
        op.compute_rank(**kwargs)
        if info_tracker:
            self.update_tracker(op)

    def set_tracker_parameters(self):
        """Initialize the info tracker by setting its parameters.

        Set the names of the variables to be displayed.
        """
        self.info_tracker.set_header_list(self._get_info_header_list())

    def start_tracker(self):
        """Start the info tracker.

        Track information about the properties of the underlying operator matrices and display it in a web page.
        """
        self.info_tracker = DisplayInfo.InfoTracker(str(self))
        self.set_tracker_parameters()
        op_info_dict = collections.OrderedDict()
        for op in self.op_matrix_list:
            op_info_dict.update({tuple(
                op.get_domain().get_ordered_param_dict().values()): op.get_properties().list()})
        self.info_tracker.update_data(op_info_dict)
        self.info_tracker.start()

    def update_tracker(self, op):
        """Update info tracker for the operator matrix op.

        :param op: Operator matrix for which to update the properties and message it to the info tracker.
        :type op: OperatorMatrix
        """
        op.update_properties()
        message = {tuple(op.get_domain().get_ordered_param_dict(
        ).values()): op.get_properties().list()}
        self.info_tracker.get_queue().put(message)

    def stop_tracker(self):
        """Stop tracking information about the underlying operator matrices."""
        self.info_tracker.stop()

    def plot_info(self):
        """Plot information about the operator matrices of the operator matrix collection to a html file."""
        path = self.get_info_plot_path()
        header_list = self._get_info_header_list()
        data_list = []
        for op in self.op_matrix_list:
            op.update_properties()
            data_list.append(
                list(op.domain.get_ordered_param_dict().values()) + op.get_properties().list())
        DisplayInfo.plot_info(data_list, header_list,
                              path, to_html=True, to_csv=False)

    def _get_info_header_list(self):
        try:
            param_names = list(self.get_vector_space().get_vs_list()[
                               0].get_ordered_param_dict().keys())
        except IndexError:
            param_names = []
        return param_names + OperatorMatrixProperties.names()


class Differential(OperatorMatrixCollection):
    """Operator matrix collection, which is supposed to be a differential, i.e. is supposed to square to zero."""
    __metaclass__ = ABCMeta

    def __init__(self, sum_vector_space, op_matrix_list):
        """Initialze the underlying sum vector space and the list of operator matrices composing the differential.

        :param sum_vector_space: Underlying vector space.
        :type sum_vector_space: GraphVectorSpace.SumVectorSpace
        :param op_matrix_list: List of operator matrices composing the differential
        :type op_matrix_list: list(OperatorMatrix)
        """
        super(Differential, self).__init__(sum_vector_space, op_matrix_list)

    @abstractmethod
    def get_cohomology_plot_path(self):
        """Return the path to the cohomology plot file.

        File name without ending.

        :return: Path to the cohomology plot file.
        :rtype: path
        """
        pass

    @abstractmethod
    def get_cohomology_web_path(self):
        """Return the path to the data file for the web.

        File name without ending.

        :return: Path to the data file.
        :rtype: path
        """
        pass

    def get_cohomology_plot_parameter_order(self):
        pass

    def get_ordered_cohomology_param_range_dict(self):
        """Return an ordered dictionary of parameter ranges for the sub vector spaces of the underlying sum vector space.

         :return: Ordered dictionary of parameter ranges.
         :rtype: Shared.OrderedDict
         :Example: Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range)])
         """
        return self.sum_vector_space.get_ordered_param_range_dict()

    def __str__(self):
        """Return a unique description of the differential.

        :return: Unique description of the differential.
        :rtype: str
        """
        return '<%s differential on %s>' % (self.get_type(), str(self.sum_vector_space))

    @staticmethod
    def cohomology_dim(opD, opDD):
        """Compute the cohomology dimension, i.e. dim(ker(opD)/im(opDD)) = dim(opD.domain) - rank(opD) - rank(opDD).

        The domain vector space of opD must correspond to the target vector space of opDD.
        A warning is printed a negative cohomology dimension results.

        :param opD: First operator matrix relevant for the cohomology dimension.
        :type opD: OperatorMatrix
        :param opDD: Second operator matrix relevant for the cohomology dimension.
        :type opDD: OperatorMatrix
        :return: Cohomology dimension dim(ker(opD)/im(opDD)) = dim(opD.domain) - rank(opD) - rank(opDD).
        :rtype: int
        """
        try:
            dimV = opD.get_domain().get_dimension()
        except StoreLoad.FileNotFoundError:
            logger.info(
                "Cannot compute cohomology: First build basis for %s " % str(opD.get_domain()))
            return None
        if dimV == 0:
            return '*'
        if opD.is_valid():
            try:
                rankD = opD.get_matrix_rank()
            except StoreLoad.FileNotFoundError:
                logger.info(
                    "Cannot compute cohomology: Matrix rank not calculated for %s " % str(opD))
                return None
        else:
            rankD = 0
        if opDD.is_valid():
            try:
                rankDD = opDD.get_matrix_rank()
            except StoreLoad.FileNotFoundError:
                logger.info(
                    "Cannot compute cohomology: Matrix rank not calculated for %s " % str(opDD))
                return None
        else:
            rankDD = 0
        cohomology_dim = dimV - rankD - rankDD
        if cohomology_dim < 0:
            raise ValueError("Negative cohomology dimension for %s (%d - %d - %d)" %
                             (str(opD.domain), dimV, rankD, rankDD))
            #logger.error("Negative cohomology dimension for %s" % str(opD.domain))
        return cohomology_dim

    def _get_cohomology_dim_dict(self):
        """Return a dictionary for the cohomology dimensions of the sub vector spaces of the underlying vector space.

        :return: Dictionary (vector space -> cohomology dimension)
        :rtype: dict(GraphVectorSpace.VectorSpace -> int)
        """
        cohomology_dim_dict = dict()
        for (opD, opDD) in itertools.permutations(self.op_matrix_list, 2):
            if opD.get_domain() == opDD.get_target():
                dim = Differential.cohomology_dim(opD, opDD)
                cohomology_dim_dict.update({opD.domain: dim})
        return cohomology_dim_dict

    def get_cohomology_dim_dict(self):
        """Return a dictionary for the cohomology dimensions of the sub vector spaces of the underlying vector space.

        :return: Dictionary (vector space parameters -> cohomology dimension)
        :rtype: dict(tuple(int) -> int)
        """
        cohomology_dim = self._get_cohomology_dim_dict()
        dim_dict = dict()
        for vs in self.sum_vector_space.get_vs_list():
            dim_dict.update(
                {vs.get_ordered_param_dict().get_value_tuple(): cohomology_dim.get(vs)})
        return dim_dict

    def complex_is_acyclic(self):
        """Test whether the associated comples is acyclic.

        :return: True if the associated comples is acyclic.
        :rtype: bool
        """
        cohomology_dims = self._get_cohomology_dim_dict().values()
        for dim in cohomology_dims:
            if dim is not None and not (dim == 0 or dim == '*'):
                return False
        return True

    def square_zero_test(self, eps=Parameters.square_zero_test_eps):
        """Generic test whether the differential squares to zero.

        Searche for matching pairs in the list of underlying operator matrices and test whether they square to zero.
        Report for how many of them the test was trivially successful (because at least two matrices are trivial),
        successful, inconclusive (because matrices are missing) or unsuccessful.

        :param eps: Threshold for the square zero test (Default: Parameters.square_zero_test_eps).
        :type eps: float
        :return: Tuple with the number of pairs for which the square zero test was a
            (trivial success, success, inconclusive, fail).
        :rtype: tuple(int, int, int, int)
        """
        print(' ')
        print("Square zero test for %s:" % str(self))
        succ_l = 0  # number of pairs for which test was successful
        fail = []  # failed pairs
        triv_l = 0  # number of pairs for which test trivially succeeded because at least one operator is the empty matrix
        inc_l = 0  # number of pairs for which operator matrices are missing
        for (op1, op2) in itertools.permutations(self.op_matrix_list, 2):
            if op1.get_target() == op2.get_domain():
                # A composable pair is found

                pair = (op1, op2)
                res = self._square_zero_test_for_pair(pair, eps=eps)

                if res == 'triv':
                    triv_l += 1
                elif res == 'succ':
                    print('success')
                    succ_l += 1
                elif res == 'inc':
                    inc_l += 1
                elif res == 'fail':
                    print("Square zero test for %s: failed for the pair %s, %s" % (
                        str(self), str(op1), str(op2)))
                    logger.error("Square zero test for %s: failed for the pair %s, %s" % (
                        str(self), str(op1), str(op2)))
                    fail.append(pair)
                else:
                    raise ValueError('Undefined commutativity test result')

        fail_l = len(fail)
        print("trivial success: %d, success: %d, inconclusive: %d, failed: %d pairs" % (
            triv_l, succ_l, inc_l, fail_l))
        logger.warninging("Square zero test for %s:" % str(self))
        logger.warninging("trivial success: %d, success: %d, inconclusive: %d, failed: %d pairs" %
                          (triv_l, succ_l, inc_l, fail_l))
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
        except StoreLoad.FileNotFoundError:
            logger.info("Cannot test square zero: "
                        "Operator matrix not built for %s or %s" % (str(op1), str(op2)))
            return 'inc'

        if Shared.matrix_norm(M2 * M1) < eps:
            return 'succ'
        return 'fail'

    def plot_cohomology_dim(self, to_html=False, to_csv=False, x_plots=2):
        """Plot the cohomology dimensions.

        Plot the cohomology dimensions as plot and/or table associated with the differential.

        :param to_html: Option to generate a html file with a table of the cohomology dimensions (Dafault: False).
        :type to_html: bool
        :param to_csv: Option to generate a csv file with a table of the cohomology dimensions (default: False).
        :type to_csv: bool
        :param x_plots: Number of plots on the x-axis (Default: 2).
        :type x_plots: int
        """
        print(' ')
        print('Plot cohomology dimensions of the associated graph complex of ' + str(self))
        logger.warning(
            'Plot cohomology dimensions of the associated graph complex of ' + str(self))
        dim_dict = self.get_cohomology_dim_dict()
        plot_path = self.get_cohomology_plot_path()
        parameter_order = self.get_cohomology_plot_parameter_order()
        ordered_param_range_dict = self.get_ordered_cohomology_param_range_dict()
        PlotCohomology.plot_array(dim_dict, ordered_param_range_dict, plot_path, to_html=to_html, to_csv=to_csv,
                                  x_plots=x_plots, parameter_order=parameter_order)

    def export_cohomology_dim_for_web(self):
        """Writes the cohomology data into a js file for use on the website."""
        export_path = self.get_cohomology_web_path() + ".js"
        StoreLoad.generate_path(export_path)
        print(' ')
        print(
            'Exporting cohomology dimensions of the associated graph complex of ' + str(self) + ' to file '+export_path+'.')
        dim_dict = self.get_cohomology_dim_dict()
        s = "data_sources.push({ name : 'GH', description : 'GH computation', data: ["
        for k, v in dim_dict.items():
            if v is not None:
                s += str(list(k) + [v]) + ",\n"
        s += "] });"
        # s = s.replace('None', 'null')

        with open(export_path, 'w') as outfile:
            outfile.write(s)
