import scipy.sparse as sparse
import logging
from sage.all import *
import Shared as SH

reload(SH)


def exists_file(path):
    return (path is not None) and os.path.isfile(path)


def _load_basis_g6(path):
    logging.info("Load basis from reference file: %s" % str(path))
    basis_g6 = SH.load_string_list(path, header=False)
    return basis_g6


def g6_to_canon_g6(graph6):
    graph = Graph(graph6)
    return graph.canonical_label().graph6_string()


def get_basis_g6(path):
    logging.info("Load basis from reference file: %s" % str(path))
    if not exists_file(path):
        raise SH.RefError("%s: Reference basis file not found" % str(path))
    basis_g6 = _load_basis_g6(path)
    basis_g6_canon = []
    for G6 in basis_g6:
        canonG6 = g6_to_canon_g6(G6)
        basis_g6_canon.append(canonG6)
    return basis_g6_canon


def _load_matrix(path):
    logging.info("Load operator matrix from reference file: %s" % str(path))
    shape = None
    matrixList = SH.load_string_list(path, header=False)
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
                    raise ValueError("%s: Found negative matrix indices: %d %d" % (str(path), min(row), min(column)))
                if max(row) >= m or max(column) >= n:
                    raise ValueError("Matrix read from reference file %s is wrong: Index outside matrix size" % str(path))
    if shape is None:
        raise ValueError("%s: Schape of reference matrix is unknown" % str(path))
    return ((data, (row, column)), shape)


def get_matrix(path):
    if not exists_file(path):
        raise SH.RefError("%s: Reference basis file not found" % str(path))
    (matrixList, shape) = _load_matrix(path)
    (m, n) = shape
    logging.info("Get reference operator matrix from file %s with shape (%d, %d)" % (str(path), m, n))
    return sparse.csr_matrix(matrixList, shape=shape, dtype=float)