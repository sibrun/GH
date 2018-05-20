import logging
from sage.all import *
import StoreLoad


class RefGraphVectorSpace(object):
    def __init__(self, graph_vs):
        self.graph_vs = graph_vs
        self.basis_file_path = graph_vs.get_ref_basis_file_path()

    def __str__(self):
        return "<Reference vector space: %s>" % str(self.basis_file_path)

    def __eq__(self, other):
        return self.graph_vs == other.graph_vs

    def exists_basis_file(self):
        return (self.basis_file_path is not None) and os.path.isfile(self.basis_file_path)

    def _load_basis_g6(self):
        if not self.exists_basis_file():
            raise StoreLoad.FileNotFoundError("%s: Reference basis file not found" % str(self))
        return StoreLoad.load_string_list(self.basis_file_path)

    def _g6_to_canon_g6(self, graph6, sign=False):
        graph = Graph(graph6)
        if not sign:
            return graph.canonical_label().graph6_string()
        canonG, perm_dict = graph.canonical_label(partition=self.graph_vs.get_partition(), certificate=True)
        sgn = self.graph_vs.perm_sign(graph, perm_dict.values())
        return (canonG.graph6_string(), sgn)

    def get_basis_g6(self):
        basis_g6 = self._load_basis_g6()
        basis_g6_canon = []
        for G6 in basis_g6:
            basis_g6_canon.append(self._g6_to_canon_g6(G6))
        return basis_g6_canon

    def get_dimension(self):
        try:
            dim = len(self._load_basis_g6())
        except StoreLoad.FileNotFoundError:
            if not self.graph_vs.is_valid():
                return 0
            else:
                raise StoreLoad.FileNotFoundError('Dimension of reference basis unknown: %s' % str(self))
        return dim

    def get_transformation_matrix(self):
        basis_g6 = self._load_basis_g6()
        dim = len(basis_g6)
        if not dim == self.graph_vs.get_dimension():
            raise ValueError("Dimension of reference basis and basis not equal for %s" % str(self))
        T = matrix(ZZ, dim, dim, sparse=True)
        lookup = {G6: j for (j, G6) in enumerate(self.graph_vs.get_basis_g6())}
        i = 0
        for G6 in basis_g6:
            (canonG6, sgn) = self._g6_to_canon_g6(G6, sign=True)
            j = lookup.get(canonG6)
            if j is None:
                raise ValueError("%s: Graph from ref basis not found in basis" % str(self))
            T.add_to_entry(i, j, sgn)
            i += 1
        if not T.is_invertible():
            raise ValueError("%s: Basis transformation matrix not invertible" % str(self))
        return T


class RefOperatorMatrix(object):
    def __init__(self, op_matrix):
        self.op_matrix = op_matrix
        self.domain = RefGraphVectorSpace(self.op_matrix.get_domain())
        self.target = RefGraphVectorSpace(self.op_matrix.get_target())
        self.matrix_file_path = self.op_matrix.get_ref_matrix_file_path()
        self.rank_file_path = self.op_matrix.get_ref_rank_file_path()

    def __str__(self):
        return "<Reference operator: %s>" % str(self.matrix_file_path)

    def exists_matrix_file(self):
        return os.path.isfile(self.matrix_file_path)

    def _load_matrix(self):
        if not self.exists_matrix_file():
           raise StoreLoad.FileNotFoundError("%s: Reference basis file not found" % str(self))
        logging.info("Load operator matrix from reference file: %s" % str(self.matrix_file_path))
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
        (entriesList, shape) = self._load_matrix()
        if shape is None:
            (d, t) = (self.domain.get_dimension(), self.target.get_dimension())
        else:
            (d, t) = shape
        logging.info("Get reference operator matrix from file %s" % str(self))
        M = matrix(ZZ, d, t, sparse=True)
        for (i, j, v) in entriesList:
            M.add_to_entry(i, j, v)
        return M.transpose()

    def get_matrix(self):
        M = self.get_matrix_wrt_ref()
        T_domain = self.domain.get_transformation_matrix()
        T_target = self.target.get_transformation_matrix()
        (m, n) = (M.nrows(), M.ncols())
        if m == 0 or n == 0:
            return M
        return T_target.inverse() * M * T_domain

    def exists_rank_file(self):
        return os.path.isfile(self.rank_file_path)

    def get_rank(self):
        if not self.exists_rank_file():
           raise StoreLoad.FileNotFoundError("%s: Reference rank file not found" % str(self))
        return int(StoreLoad.load_line(self.rank_file_path))
