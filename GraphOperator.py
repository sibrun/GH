from abc import ABCMeta, abstractmethod
import os
import scipy.sparse as sparse
import numpy as np
import pickle
import logging
from sage.all import *
import Shared as SH

reload(SH)

class GraphOperator():
    __metaclass__ = ABCMeta

    def __init__(self, domain, target, header_ref=False, skip_if_no_basis=True):
        self.domain = domain
        self.target = target
        self.valid = self.domain.valid and self.target.valid
        self.file_path = self._set_file_path()
        self.work_estimate = self._set_work_estimate()
        self.file_path_ref = self._set_file_path(ref=True)
        self.header_ref = header_ref
        self.skip_if_no_basis = skip_if_no_basis
        self.shape = None
        self.entries = None

    @abstractmethod
    def _set_file_path(self, ref=False):
        pass

    @abstractmethod
    def _set_work_estimate(self):
        pass

    @abstractmethod
    def _operate_on(self,graph):
        pass

    @abstractmethod
    def __str__(self):
        pass

    def get_info(self):
        validity = "valid" if self.valid else "not valid"
        built = "matrix built" if self.matrix_built() else "matrix not built"
        shape = "shape unknown"
        trivial = "triviality unknown"
        entries = "entries unknown"
        if self.matrix_built():
            shape = "shape = none" if self.shape is None else "shape = (%d,%d)" % (self.shape[0],self.shape[1])
            if self.shape is not None:
                trivial = "trivial" if self.shape[0]==0 or self.shape[1]==0 else "not trivial"
            entries = "entries = none" if self.entries is None else "%d entries" % self.entries
        return "%s, %s, %s, %s, %s" % (validity, built, shape, trivial, entries)

    def build_matrix(self):
        if not self.valid:
            logging.info("Skip creating file for operator matrix, since %s is not valid" % str(self))
            return
        if self.matrix_built():
            self.load_info_from_file()
            return
        try:
            domainBasis = self.domain.get_basis(g6=False)
        except SH.NotBuiltError:
            if not self.skip_if_no_basis:
                raise SH.NotBuiltError("Cannot build operator matrix of %s: First build basis of the domain %s" % (str(self), str(self.domain)))
            else:
                logging.warn("Skip building operator matrix of %s since basis of the domain %s is not built" % (str(self), str(self.domain)))
                return
        try:
            targetBasis6 = self.target.get_basis(g6=True)
        except SH.NotBuiltError:
            if not self.skip_if_no_basis:
                raise SH.NotBuiltError("Cannot build operator matrix of %s: First build basis of the target %s" % (str(self), str(self.target)))
            else:
                logging.warn("Skip building operator matrix of %s since basis of the target %s is not built" % (str(self), str(self.target)))
                return

        self.shape = (self.domain.dimension, self.target.dimension)
        if self.is_trivial():
            self.entries = 0
            self._store_matrix([])
            logging.info("Created file with empty operator matrix of %s, since it's shape is %s" % (str(self), str(self.shape)))
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
                imageIndex = lookup.get(image)
                if imageIndex is not None:
                    matrixList.append("%d %d %d" % (domainIndex, imageIndex, factor))
        self.entries = len(matrixList)
        self._store_matrix(matrixList)
        logging.info("Operator matrix built for %s" % str(self))

    def matrix_built(self):
        if os.path.isfile(self.file_path):
            return True
        return False

    def load_info_from_file(self):
        try:
            header = SH.load_header(self.file_path)
        except SH.FileNotExistingError:
            logging.warn("Cannot load infos from file %s" % str(self.file_path))
            return
        self._load_info_from_header(header)
        self.domain.load_info_from_file()
        self.target.load_info_from_file()
        logging.info("Infos loaded from file %s" % str(self.file_path))

    def _load_info_from_header(self, header):
        (m, n, self.entries) = map(int,header.split(" "))
        self.shape = (m, n)

    def _get_header(self):
        (m, n) = self.shape
        return "%d %d %d" % (m, n, self.entries)

    def is_trivial(self):
        if self.shape is None:
            raise SH.NotBuiltError("Shape of %s unknown: First load matrix" % str(self))
        return self.shape[0] == 0 or self.shape[1] == 0

    def _store_matrix(self, matrixList):
        logging.info("Store operator matrix in file: %s" % str(self.file_path))
        SH.store_string_list(matrixList, self.file_path, header=self._get_header())

    def _load_matrix(self):
        logging.info("Load operator matrix from file: %s" % str(self.file_path))
        (header, matrixList) = SH.load_string_list(self.file_path, header=True)
        self._load_info_from_header(header)
        if len(matrixList) != self.entries:
            raise ValueError("Number of matrix entries read from file %s is wrong" % str(self.file_path))
        row = []
        column = []
        data = []
        for line in matrixList:
            (i, j, v) = map(int, line.split(" "))
            row.append(i)
            column.append(j)
            data.append(v)
        return (row, column, data)

    def _load_ref_matrix(self):
        logging.info("Load operator matrix from reference file: %s" % str(self.file_path_ref))
        data = SH.load_string_list(self.file_path_ref, header=self.header_ref)
        if self.header_ref:
            (header, matrixList) = data
            (m, n, entries) = map(int, header.split(" "))
            if len(matrixList) != entries:
                raise ValueError("Number of matrix entries read from file %s is wrong" % str(self.file_path_ref))
        else:
            matrixList = data
        row = []
        column = []
        data = []
        for line in matrixList:
            (i, j, v) = map(int, line.split(" "))
            row.append(i-1)
            column.append(j-1)
            data.append(v)
        if data[-1] == 0:
            if max(row) > row[-1] or max(column) > column[-1]:
                raise ValueError("Matrix read from referencde file %s is wrong: Index exceeds matrix diemension" % str(self.file_path_ref))
        return (row, column, data)

    def get_matrix(self):
        if not self.valid:
            log.warning("No operator matrix: %s is not valid" % str(self))
            return None
        if not self.matrix_built():
            raise SH.NotBuiltError("Cannot load operator matrix, No operator file found for %s: " % str(self))
        matrix = self._load_matrix()
        logging.info("Get operator matrix of %s with shape %s" % (str(self), str(self.shape)))
        (row, column, data) = matrix
        return sparse.csr_matrix((data, (row, column)), shape=self.shape, dtype=float)

    def ref_file_available(self):
        return (self.file_path_ref is not None) and os.path.isfile(self.file_path_ref)

    def get_ref_matrix(self):
        if self.file_path_ref is None:
            raise SH.RefError("%s: Path to reference file not specified" % str(self))
        if not os.path.isfile(self.file_path_ref):
            raise SH.RefError("%s: Reference basis file not found" % str(self))
        (row, column, data) = self._load_ref_matrix()
        if len(row):
            shape=(row.pop(-1)+1,column.pop(-1)+1)
            data.pop(-1)
        else:
            shape=(0,0)
        return sparse.csr_matrix((data, (row, column)), shape=shape, dtype=float)

    def delete_file(self):
        if os.path.isfile(self.file_path):
            os.remove(self.file_path)
