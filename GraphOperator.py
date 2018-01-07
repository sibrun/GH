from abc import ABCMeta, abstractmethod
import itertools
import scipy.sparse as sparse
import logging
from sage.all import *
import Shared as SH

reload(SH)


class GraphOperator():
    __metaclass__ = ABCMeta

    def __init__(self, domain, target):
        self.domain = domain
        self.target = target
        self.valid = self.domain.valid and self.target.valid
        self.file_path = self._set_file_path()

    @classmethod
    @abstractmethod
    def generate_operators(cls, vs_list):
        pass

    @abstractmethod
    def _set_file_path(self):
        pass

    @abstractmethod
    def get_file_path_ref(self):
        pass

    @abstractmethod
    def get_work_estimate(self):
        pass

    @abstractmethod
    def _operate_on(self,graph):
        pass

    @abstractmethod
    def __str__(self):
        pass

    def get_info(self):
        validity = "valid" if self.valid else "not valid"
        built = "matrix built" if self.exists_file() else "matrix not built"
        shape = "matrix shape unknown"
        entries = "entries unknown"
        if self.exists_file():
            (s, e) = self.get_matrix_header()
            (m, n) = s
            shape = "matrix shape = (%d, %d)" % (m, n)
            entries = "%d entries" % e
        return "%s, %s, %s" % (validity, shape, entries)

    def build_matrix(self, skip_if_no_basis = True):
        if not self.valid:
            logging.info("Skip creating file for operator matrix, since %s is not valid" % str(self))
            return
        if self.exists_file():
            return
        try:
            domainBasis = self.domain.get_basis(g6=False)
        except SH.NotBuiltError:
            if not skip_if_no_basis:
                raise SH.NotBuiltError("Cannot build operator matrix of %s: First build basis of the domain %s" % (str(self), str(self.domain)))
            else:
                logging.warn("Skip building operator matrix of %s since basis of the domain %s is not built" % (str(self), str(self.domain)))
                return
        try:
            targetBasis6 = self.target.get_basis(g6=True)
        except SH.NotBuiltError:
            if not skip_if_no_basis:
                raise SH.NotBuiltError("Cannot build operator matrix of %s: First build basis of the target %s" % (str(self), str(self.target)))
            else:
                logging.warn("Skip building operator matrix of %s since basis of the target %s is not built" % (str(self), str(self.target)))
                return

        shape = (self.target.get_dimension(), self.domain.get_dimension())
        (m, n) = shape
        if m == 0 or n == 0:
            entries = 0
            self._store_matrix([],(shape, entries))
            logging.info("Created empty file for operator matrix of %s, since the matrix shape is (%d, %d)" % (str(self), m, n))
            return

        lookup = {s: j for (j,s) in enumerate(targetBasis6)}
        matrixList = []
        for (domainIndex, G) in enumerate(domainBasis):
            imageList = self._operate_on(G)
            canonImages = dict()
            for (GG, prefactor) in imageList:
                (GGcanon6, sgn1) = self.target.graph_to_canon_g6(GG)
                sgn0 = canonImages.get(GGcanon6)
                if sgn0 is None:
                    sgn0 = 0
                canonImages.update({GGcanon6: (sgn0 + sgn1 * prefactor)})
            for (image, factor) in canonImages.items():
                if factor:
                    targetIndex = lookup.get(image)
                    if targetIndex is not None:
                        matrixList.append("%d %d %d" % (targetIndex, domainIndex, factor))
        entries = len(matrixList)
        self._store_matrix(matrixList,(shape, entries))
        logging.info("Operator matrix built for %s" % str(self))

    def exists_file(self):
        return os.path.isfile(self.file_path)

    def exist_domain_target_files(self):
        return os.path.isfile(self.domain.file_path) and os.path.isfile(self.target.file_path)

    def get_matrix_header(self):
        if not self.valid:
            return ((self.target.get_dimension(), self.domain.get_dimension()), 0)
        try:
            header = SH.load_header(self.file_path)
        except SH.FileNotExistingError:
            raise SH.NotBuiltError("Cannot load header from file %s: Build operator matrix first" % str(self.file_path))
        return GraphOperator._get_matrix_info_from_header(header)

    def get_matrix_entries(self):
        (shape, entries) = self.get_matrix_header()
        return entries

    @staticmethod
    def _get_matrix_info_from_header(header):
        (m, n, entries) = map(int, header.split(" "))
        shape = (m, n)
        return(shape, entries)

    @staticmethod
    def _get_header_from_matrix_info(matrix_info):
        (shape, entries) = matrix_info
        (m, n) = shape
        return "%d %d %d" % (m, n, entries)

    def is_trivial(self):
        if not self.valid:
            return True
        try:
            (shape, entries) = self.get_matrix_header()
            (m, n) = shape
        except SH.NotBuiltError:
            try:
                m = self.domain.get_dimension()
                n = self.target.get_dimension()
                entries = 1
            except SH.NotBuiltError:
                raise SH.NotBuiltError("Matrix shape of %s unknown: Build operator matrix first" % str(self))
        return m == 0 or n == 0 or entries == 0

    def _store_matrix(self, matrixList, matrix_info):
        logging.info("Store operator matrix in file: %s" % str(self.file_path))
        header = GraphOperator._get_header_from_matrix_info(matrix_info)
        SH.store_string_list(matrixList, self.file_path, header=header)

    def _load_matrix(self):
        logging.info("Load operator matrix from file: %s" % str(self.file_path))
        (header, matrixList) = SH.load_string_list(self.file_path, header=True)
        (shape, entries) = GraphOperator._get_matrix_info_from_header(header)
        (m, n) = shape
        if m != self.target.get_dimension() or n != self.domain.get_dimension():
            raise ValueError("Shape of matrix doesn't correspond to the vector space dimensions for the file: %s" % str(self.file_path))
        if len(matrixList) != entries:
            raise ValueError("Number of matrix entries read from file %s is wrong" % str(self.file_path))
        row = []
        column = []
        data = []
        for line in matrixList:
            (i, j, v) = map(int, line.split(" "))
            if not (i >= 0 and j >= 0):
                raise ValueError("%s: Found negative matrix indices: %d %d" % (str(self.file_path),i, j))
            row.append(i)
            column.append(j)
            data.append(v)
        if len(row):
            if min(row) < 0 or min(column) < 0:
                raise ValueError("%s: Found negative matrix indices: %d %d" % (str(self.file_path), min(row), min(column)))
            if max(row) >= m or max(column) >= n:
                raise ValueError("Matrix read from file %s is wrong: Index outside matrix size" % str(self.file_path))
        return ((data, (row, column)), shape)

    def get_matrix(self):
        if not self.valid:
            log.warning("No operator matrix: %s is not valid" % str(self))
            shape = (m ,n) = (self.target.get_dimension(), self.domain.get_dimension())
            matrixList = ([], ([], []))
        else:
            if not self.exists_file():
                raise SH.NotBuiltError("Cannot load operator matrix, No operator file found for %s: " % str(self))
            (matrixList, shape) = self._load_matrix()
            (m, n) = shape
        logging.info("Get operator matrix of %s with shape (%d, %d)" % (str(self), m, n))
        matrix = sparse.csc_matrix(matrixList, shape=shape, dtype=int)
        matrix.eliminate_zeros()
        return matrix

    def delete_file(self):
        if os.path.isfile(self.file_path):
            os.remove(self.file_path)



