"""Module providing implementations for graph vector spaces and operator matrices based on reference data.

It provides functionality to compare the basis of a vector space as well as the operator matrices with
reference data.
"""

__all__ = ['RefGraphVectorSpace', 'RefOperatorMatrix']

import Log
from sage.all import *
import StoreLoad


logger = Log.logger.getChild('ref_graph_complex')


class RefGraphVectorSpace(object):
    """Reference graph vector space.

    Reads the graph vector space basis from the reference data file and provides
    the transformation matrix to change the coordinate from the reference basis to the basis.

    Attributes:
        graph_vs (GraphVectorSpace.GraphVectorSpace): Graph vector space to be compared with reference data.

        basis_file_path (path): Path to reference basis file.
    """
    def __init__(self, graph_vs):
        """Initialize graph vector space.

        :param graph_vs: GraphVectorSpace.GraphVectorSpace: Graph vector space to be compared with reference data.
        """
        self.graph_vs = graph_vs
        self.basis_file_path = graph_vs.get_ref_basis_file_path()

    def __str__(self):
        """Returns a unique description of the reference graph vector space based on the reference file path.

        :return: str: Description of the reference graph vector space.
        """
        return "<Reference vector space: %s>" % str(self.basis_file_path)

    def __eq__(self, other):
        """Compares two reference graph vector spaces.

        :param other: RefGraphVectorSpace: Reference graph vector space to be compared with.
        :return: bool: True if the two reference graph vector spaces are equal.
        """
        return self.graph_vs == other.graph_vs

    def exists_basis_file(self):
        """Checks whether the reference basis file exists.

        :return: bool: True if the reference basis file is found.
        """
        return (self.basis_file_path is not None) and os.path.isfile(self.basis_file_path)

    def _load_basis_g6(self):
        """ Loads the reference basis list from the reference file.

        The implementation depends on the reference data.

        :return: list(str): List of graph6 strings representing the reference basis.
        :raise StoreLoad.FileNotFoundError: If the reference basis file is not found.
        """
        if not self.exists_basis_file():
            raise StoreLoad.FileNotFoundError("%s: Reference basis file not found" % str(self))
        return StoreLoad.load_string_list(self.basis_file_path)

    def _g6_to_canon_g6(self, graph6, sgn=False):
        graph = Graph(graph6)
        if not sgn:
            return graph.canonical_label(partition=self.graph_vs.get_partition()).graph6_string()
        canonG, perm_dict = graph.canonical_label(partition=self.graph_vs.get_partition(), certificate=True)
        sgn = self.graph_vs.perm_sign(graph, perm_dict.values())
        return (canonG.graph6_string(), sgn)

    def get_basis_g6(self):
        """Returns the reference basis as list of graph6 strings.

        :return: list(str): Reference basis as list of graph6 strings.
        :raise StoreLoad.FileNotFoundError: If the reference basis file is not found.
        """
        basis_g6 = self._load_basis_g6()
        basis_g6_canon = []
        for G6 in basis_g6:
            basis_g6_canon.append(self._g6_to_canon_g6(G6))
        return basis_g6_canon

    def get_dimension(self):
        """Returns the dimension of the reference vector space.

        :return: non-negative int: Dimension of the reference vector space.
        :raise StoreLoad.FileNotFoundError: If the reference basis file is not found.
        """
        return len(self._load_basis_g6())

    def get_transformation_matrix(self):
        """Returns the coordinate transformation matrix, to transform from the reference basis to the basis.

        :return:
        """
        basis_g6 = self._load_basis_g6()
        dim = len(basis_g6)
        if not dim == self.graph_vs.get_dimension():
            raise ValueError("Dimension of reference basis and basis not equal for %s" % str(self))
        T = matrix(ZZ, dim, dim, sparse=True)
        lookup = {G6: j for (j, G6) in enumerate(self.graph_vs.get_basis_g6())}
        for (i, G6) in enumerate(basis_g6):
            (canonG6, sgn) = self._g6_to_canon_g6(G6, sgn=True)
            j = lookup.get(canonG6)
            if j is None:
                raise ValueError("%s: Graph from ref basis not found in basis" % str(self))
            T.add_to_entry(i, j, sgn)
        if not T.is_invertible():
            raise ValueError("%s: Basis transformation matrix not invertible" % str(self))
        return T


class RefOperatorMatrix(object):
    """Reference operator matrix.

    Reads the operator matrix from the reference data file and provides operator matrix w.r.t. the reference basis
    as well as to the basis.

    Attributes:
        op_matrix (GraphOperator.GraphOperatorMatrix): Graph operator matrix to be compared with reference data.

        domain (RefGraphVectorSpace): Reference domain vector space of the reference matrix.

        target (RefGraphVectorSpace): Reference target vector space of the reference matrix.

        matrix_file_path (path): Path to the reference matrix file.

        rank_file_path (path): Path to the reference rank file.

    """
    def __init__(self, op_matrix):
        self.op_matrix = op_matrix
        self.domain = RefGraphVectorSpace(self.op_matrix.get_domain())
        self.target = RefGraphVectorSpace(self.op_matrix.get_target())
        self.matrix_file_path = self.op_matrix.get_ref_matrix_file_path()
        self.rank_file_path = self.op_matrix.get_ref_rank_file_path()

    def __str__(self):
        """Returns a unique description of the reference operator matrix based on the reference file path.

        :return: str: Description of the reference operator matrix.
        """
        return "<Reference operator: %s>" % str(self.matrix_file_path)

    def exists_matrix_file(self):
        """Checks whether the reference matrix file exists.

        :return: bool: True if the reference matrix file is found.
        """
        return os.path.isfile(self.matrix_file_path)

    def _load_matrix(self):
        """ Loads the reference matrix list from the reference file.

        The implementation depends on the reference data.

        :return: tuple(list(tuple(non-negative int, non-negative int, int)), tuple(non-negative int, non-negative int)):
            (matrix_list = list((domain index, target index, value), shape = (domain dimension, target dimension))
        :raise StoreLoad.FileNotFoundError: If the reference matrix file is not found.
        :raise: ValueError: Raised in the following cases:
            End line missing.
            Negative matrix index.
            Matrix index too large.
        """
        if not self.exists_matrix_file():
           raise StoreLoad.FileNotFoundError("%s: Reference basis file not found" % str(self))
        stringList = StoreLoad.load_string_list(self.matrix_file_path)
        entriesList = []
        if len(stringList) == 0:
            return ([], None)
        else:
            (d, t, z) = map(int, stringList.pop().split(" "))
            if z != 0:
                raise ValueError("End line in reference file %s is missing" % str(self))
            shape = (d, t)
            for line in stringList:
                (i, j, v) = map(int, line.split(" "))
                if i < 0 or j < 0:
                    raise ValueError("%s: Negative matrix index" % str(self))
                if i > d or j > t:
                    raise ValueError("%s Matrix index outside matrix size" % str(self))
                if i == 0 or j == 0:
                    continue
                entriesList.append((i - 1, j - 1, v))
        return (entriesList, shape)

    def get_matrix_wrt_ref(self):
        """Returns the reference matrix w.r.t. the reference basis.

        :return: sage.Matrix_sparse: Reference Matrix as sparse sage matrix w.r.t. the reference basis.
        :raise StoreLoad.FileNotFoundError: If the reference matrix file or reference basis is not found.
        :raise: ValueError: Raised in the following cases:
            End line missing.
            Negative matrix index.
            Matrix index too large.
        """

        (entriesList, shape) = self._load_matrix()
        if shape is None:
            (d, t) = (self.domain.get_dimension(), self.target.get_dimension())
        else:
            (d, t) = shape
        M = matrix(ZZ, d, t, sparse=True)
        for (i, j, v) in entriesList:
            M.add_to_entry(i, j, v)
        return M.transpose()

    def get_matrix(self):
        """Returns the reference matrix w.r.t. the basis.

        :return: sage.Matrix_sparse: Reference Matrix as sparse sage matrix w.r.t. the basis.
        :raise StoreLoad.FileNotFoundError: If the reference matrix file or reference basis is not found.
        :raise: ValueError: Raised in the following cases:
            End line missing.
            Negative matrix index.
            Matrix index too large.
           """
        M = self.get_matrix_wrt_ref()
        T_domain = self.domain.get_transformation_matrix()
        T_target = self.target.get_transformation_matrix()
        (m, n) = (M.nrows(), M.ncols())
        if m == 0 or n == 0:
            return M
        return T_target.inverse() * M * T_domain

    def exists_rank_file(self):
        """Checks whether the reference rank file exists.

        :return: bool: True if the reference rank file is found.
        """
        return os.path.isfile(self.rank_file_path)

    def get_rank(self):
        """Returns the reference rank from the reference rank file.

        :return: non-negative int: Reference rank.
        :raise StoreLoad.FileNotFoundError: If refrence rank file not found.
        """
        if not self.exists_rank_file():
           raise StoreLoad.FileNotFoundError("%s: Reference rank file not found" % str(self))
        return int(StoreLoad.load_line(self.rank_file_path))
