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
    def __init__(self, vs):
        self.vs = vs
        self.file_path_ref = vs.file_path_ref

    def __str__(self):
        return "Reference vector space: %s" % str(self.file_path_ref)

    def exists_file(self):
        return (self.file_path_ref is not None) and os.path.isfile(self.file_path_ref)

    def _load_basis_g6(self):
        logging.info("Load basis from reference file: %s" % str(self))
        basis_g6 = SH.load_string_list(self.file_path_ref, header=False)
        return basis_g6

    def g6_to_canon_g6(self, graph6, sgn=False):
        graph = Graph(graph6)
        if not sgn:
            return graph.canonical_label().graph6_string()
        canonG, permDict = graph.canonical_label(certificate=True)
        sgn = self.vs.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sgn)

    def get_basis_g6(self):
        if not self.exists_file():
            raise SH.RefError("%s: Reference basis file not found" % str(self))
        basis_g6 = self._load_basis_g6()
        basis_g6_canon = []
        for G6 in basis_g6:
            basis_g6_canon.append(self.g6_to_canon_g6(G6))
            return basis_g6_canon

    def build_transformation_matrix(self):
        if not self.exists_file():
            raise SH.RefError("%s: Reference basis file not found" % str(self))
        basis_g6 = self._load_basis_g6()
        matrixList = []
        lookup = {s: j for (j, s) in enumerate(self.vs.get_basis(g6=True))}
        i = 0
        for G6 in basis_g6:
            (canonG6, sgn) = self.g6_to_canon_g6(G6, sgn=True)
            j = lookup.get(canonG6)
            if j is None:
                raise ValueError("%s: Graph from ref basis not found in basis" % str(self))
            matrixList.append("%d %d %d" % (i, j, sgn))


class RefOperator:
    def __init__(self, op):
        self.op = op
        self.file_path_ref = op.file_path_ref

    def __str__(self):
        return "Reference operator: %s" % str(self.file_path_ref)

    def exists_file(self):
        return (self.file_path_ref is not None) and os.path.isfile(self.file_path_ref)

    def _load_matrix(self):
        logging.info("Load operator matrix from reference file: %s" % str(self.file_path_ref))
        shape = None
        matrixList = SH.load_string_list(self.file_path_ref, header=False)
        row = []
        column = []
        data = []
        if matrixList is []:
            shape = (0,0)
        else:
            for line in matrixList:
                (i, j, v) = map(int, line.split(" "))
                row.append(i-1)
                column.append(j-1)
                data.append(v)
            if data[-1] == 0:
                shape = (m, n) = (row.pop(-1)+1, column.pop(-1)+1)
                data.pop(-1)
                if len(row):
                    if min(row) < 0 or min(column) < 0:
                        raise ValueError("%s: Found negative matrix indices: %d %d" % (str(self), min(row), min(column)))
                    if max(row) >= m or max(column) >= n:
                        raise ValueError("Matrix read from reference file %s is wrong: Index outside matrix size" % str(self))
        if shape is None:
            raise ValueError("%s: Shape of reference matrix is unknown" % str(self))
        return ((data, (row, column)), shape)

    def get_matrix_wrt_ref(self)
        if not self.exists_file():
           raise SH.RefError("%s: Reference basis file not found" % str(self))
        (matrixList, shape) = self._load_matrix()
        (m, n) = shape
        logging.info("Get reference operator matrix from file %s with shape (%d, %d)" % (str(self), m, n))
        return sparse.csr_matrix(matrixList, shape=shape, dtype=float)