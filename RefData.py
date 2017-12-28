from abc import ABCMeta, abstractmethod
import scipy.sparse as sparse
import logging
from sage.all import *
import Shared as SH
import GraphVectorSpace as GVS
import GraphOperator as GO

reload(SH)
reload(GVS)
reload(GO)


class RefVectorSpace:
    __metaclass__ = ABCMeta

    def __init__(self, vs):
        self.vs = vs
        self.file_path_ref = vs.get_file_path_ref()

    def __str__(self):
        return "Reference vector space: %s" % str(self.file_path_ref)

    def exists_file(self):
        return (self.file_path_ref is not None) and os.path.isfile(self.file_path_ref)

    @abstractmethod
    def _load_basis_g6(self, header):
        pass

    def _g6_to_canon_g6(self, graph6, sign=False):
        graph = Graph(graph6)
        if not sign:
            return graph.canonical_label().graph6_string()
        canonG, permDict = graph.canonical_label(certificate=True)
        sign = self.vs.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sign)

    def get_basis_g6(self, header=False):
        if not self.exists_file():
            raise SH.RefError("%s: Reference basis file not found" % str(self))
        basis_g6 = self._load_basis_g6(header)
        logging.info("Load basis from reference file: %s" % str(self))
        basis_g6_canon = []
        for G6 in basis_g6:
            basis_g6_canon.append(self._g6_to_canon_g6(G6))
            return basis_g6_canon

    def get_transformation_matrix(self, header=False):
        if not self.exists_file():
            raise SH.RefError("%s: Reference basis file not found" % str(self))
        basis_g6 = self._load_basis_g6(header)
        dim = len(basis_g6)
        if not dim == self.vs.get_dimension():
            raise ValueError("Dimension of reference basis and basis not equal for %s" % str(self))
        row = []
        column = []
        data = []
        lookup = {s: j for (j, s) in enumerate(self.vs.get_basis(g6=True))}
        i = 0
        for G6 in basis_g6:
            (canonG6, sign) = self._g6_to_canon_g6(G6, sign=True)
            j = lookup.get(canonG6)
            if j is None:
                raise ValueError("%s: Graph from ref basis not found in basis" % str(self))
            row.append(i)
            column.append(j)
            data.append(sign)
            i += 1
        return sparse.csr_matrix((data, (row, column)), shape=(dim, dim), dtype=float)


class RefOperator:
    __metaclass__ = ABCMeta

    def __init__(self, op):
        self.op = op
        self.file_path_ref = op.get_file_path_ref()

    def __str__(self):
        return "Reference operator: %s" % str(self.file_path_ref)

    def exists_file(self):
        return (self.file_path_ref is not None) and os.path.isfile(self.file_path_ref)

    @abstractmethod
    def _load_matrix(self, header):
        pass

    def get_matrix_wrt_ref(self, header=False):
        if not self.exists_file():
           raise SH.RefError("%s: Reference basis file not found" % str(self))
        (matrixList, shape) = self._load_matrix(header)
        (m, n) = shape
        logging.info("Got reference operator matrix from file %s with shape (%d, %d)" % (str(self), m, n))
        return sparse.csr_matrix(matrixList, shape=shape, dtype=float)

    def get_matrix(self, header=False):
        M = self.get_matrix_wrt_ref(header=header)
        ref_domain = RefVectorSpace(self.op.domain)
        ref_target = RefVectorSpace(self.op.target)
        T_domain = ref_domain.get_transformation_matrix(header=header)
        T_target = ref_target.get_transformation_matrix(header=header)
        return M                                                         #TODO: Implement Matrix Transformation
