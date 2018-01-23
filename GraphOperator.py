from abc import ABCMeta, abstractmethod
import itertools
import scipy.sparse as sparse
import logging
from sage.all import *
import Shared as SH

reload(SH)


class GraphOperator():
    __metaclass__ = ABCMeta

    data_type = "R"

    def __init__(self, domain, target):
        self.domain = domain
        self.target = target
        self.valid = self.domain.valid and self.target.valid
        self.file_path = self._set_file_path()
        self.rank_file_path = self._set_rank_file_path()

    @classmethod
    @abstractmethod
    def generate_operators(cls, vs_list):
        pass

    @abstractmethod
    def _set_file_path(self):
        pass

    @abstractmethod
    def _set_rank_file_path(self):
        pass

    @abstractmethod
    def get_file_path_ref(self):
        pass

    @abstractmethod
    def get_rank_file_path_ref(self):
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

        shape = (self.domain.get_dimension(), self.target.get_dimension())
        (m, n) = shape
        if m == 0 or n == 0:
            self._store_matrix([], shape)
            logging.info("Created matrix file without entries for operator matrix of %s, since the matrix shape is (%d, %d)" % (str(self), m, n))
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
                        matrixList.append((domainIndex, targetIndex, factor))
        self._store_matrix(matrixList, shape)
        logging.info("Operator matrix built for %s" % str(self))

    def exists_file(self):
        return os.path.isfile(self.file_path)

    def exists_rank_file(self):
        return os.path.isfile(self.rank_file_path)

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

    def _store_matrix(self, matrixList, shape, data_type=data_type):
        logging.info("Store operator matrix in file: %s" % str(self.file_path))
        (m, n) = shape
        stringList = []
        stringList.append("%d %d %s" % (m, n, data_type))
        for (domainIndex, targetIndex, data) in matrixList:
            stringList.append("%d %d %d" % (domainIndex + 1, targetIndex + 1, data))
        stringList.append("0 0 0")
        SH.store_string_list(stringList, self.file_path)

    def _load_matrix(self):
        logging.info("Load operator matrix from file: %s" % str(self.file_path))
        matrixList = SH.load_string_list(self.file_path)
        header = matrixList.pop(0)
        (m, n, type) = map(int, header.split(" "))
        shape = (m, n)
        if m != self.domain.get_dimension() or n != self.target.get_dimension():
            raise ValueError("Shape of matrix doesn't correspond to the vector space dimensions for the file: %s" % str(self.file_path))
        tail = map(int, matrixList.pop().split(" "))
        if not tail == (0, 0, 0):
            raise ValueError("Number of matrix entries read from file %s is wrong" % str(self.file_path))
        matrixList = []
        for line in matrixList:
            (i, j, v) = map(int, line.split(" "))
            if not (i >= 1 and j >= 1):
                raise ValueError("%s: Found invalid matrix indices: %d %d" % (str(self.file_path),i, j))
            if i >= m or j >= n:
                raise ValueError("%s: Found invalid matrix indices outside matrix size: %d %d" % (str(self.file_path),i, j))
            matrixList.append((i, j, v))
        return (matrixList, shape)

    def get_matrix(self):
        if not self.valid:
            log.warning("No operator matrix: %s is not valid" % str(self))
            shape = (m ,n) = (self.target.get_dimension(), self.domain.get_dimension())
            matrixList = ([], sahpe)
        else:
            if not self.exists_file():
                raise SH.NotBuiltError("Cannot load operator matrix, No operator file found for %s: " % str(self))
            (matrixList, shape) = self._load_matrix()
            (m, n) = shape
        logging.info("Get operator matrix of %s with shape (%d, %d)" % (str(self), m, n))
        m = matrix(ZZ, m, n, sparse=True)
        for (i, j, v) in matrixList:
            m.add_to_entry(i, j, v)
        return m

    def compute_rank(self):
        rank = self.get_matrix().rank()
        SH.store(rank, self.rank_file_path)

    def get_rank(self):
        if not self.exists_rank_file():
            raise SH.NotBuiltError("Cannot load operator rank, No rank file found for %s: " % str(self))
        SH.load(self.rank_file_path)



    def delete_file(self):
        if os.path.isfile(self.file_path):
            os.remove(self.file_path)