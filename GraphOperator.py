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
        self.file_name = self._set_file_name()
        self.work_estimate = self._set_work_estimate()
        self.file_name_ref = self._set_file_name(ref=True)
        self.header_ref = header_ref
        self.skip_if_no_basis = skip_if_no_basis
        self.shape = None
        self.entries = None

    @abstractmethod
    def _set_file_name(self, ref=False):
        pass

    @abstractmethod
    def _set_work_estimate(self):
        pass

    @abstractmethod
    def _operate_on(self,graph):
        pass

    @abstractmethod
    def params_to_string(self):
        pass

    def build_matrix(self):
        if not self.valid:
            logging.info("Skip building operator matrix for invalid parameters: " + self.params_to_stirng())
            return
        if not self.matrix_built():
            try:
                domainBasis = self.domain.get_basis(g6=False)
            except SH.NotBuiltError:
                if not self.skip_if_no_basis:
                    raise SH.NotBuiltError("Cannot build operator matrix: First build basis of the domain")
                else:
                    logging.warn("Skip building operator matrix for parameters: " + self.params_to_string() + ", since domain basis not built")
                    return
            try:
                targetBasis6 = self.target.get_basis(g6=True)
            except SH.NotBuiltError:
                if not self.skip_if_no_basis:
                    raise SH.NotBuiltError("Cannot build operator matrix: First build basis of the target")
                else:
                    logging.warn("Skip building operator matrix for parameters: " + self.params_to_string() + ", since target basis not built")
                    return

            self.shape = (self.domain.dimension, self.target.dimension)
            if self.shape[0] == 0 or self.shape[1] == 0:
                self._store_matrix([])
                logging.info("Created file with empty operator matrix for parameters: " + self.params_to_string() + ", since domain or target has dimension 0")
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
                        matrixList.append((domainIndex, imageIndex, factor))
            self._store_matrix(matrixList)
            logging.info("Created operator matrix file for parameters: " + self.params_to_string())

    def matrix_built(self):
        if os.path.isfile(self.file_name):
            return True
        return False

    def _store_matrix(self, matrixList):
        logging.info("Store operator matrix in file: " + self.file_name)
        self.entries = len(matrixList)
        out_dir = os.path.dirname(self.file_name)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        with open(self.file_name, 'w') as f:
            (m, n) = self.shape
            header = "%d %d %d" % (m, n, self.entries)
            f.write(header + '\n')
            for line in matrixList:
                (i, j, v) = line
                line_string = "%d %d %d" % (i, j, v)
                f.write(line_string + '\n')

    def _load_matrix(self):
        logging.info("Access operator file: " + self.file_name)
        with open(self.file_name, 'r') as f:
            (m, n, self.entries) = map(int,f.readline().split(" "))
            self.shape = (m, n)
            matrixList = f.read().splitlines()
        row = []
        column = []
        data = []
        for line in matrixList:
            (i, j, v) = map(int, line.split(" "))
            row.append(i)
            column.append(j)
            data.append(v)
        if len(row) != self.entries:
            raise ValueError("Number of matrix entries read from file is wrong")
        return (row, column, data)

    def _load_ref_matrix(self):
        logging.info("Access reference file: " + self.file_name_ref)
        with open(self.file_name_ref, 'r') as f:
            if self.header_ref:
                (m, n) = map(int,f.readline().split(" "))
            matrixList = f.read().splitlines()
        row = []
        column = []
        data = []
        for line in matrixList:
            (i, j, v) = map(int, line.split(" "))
            row.append(i-1)
            column.append(j-1)
            data.append(v)
        return (row, column, data)

    def get_matrix(self):
        if not self.valid:
            log.warning("No matrix for invalid parameters: " + self.params_to_string())
            return None
        if not self.matrix_built():
            raise SH.NotBuiltError("Cannot load operator matrix: No operator matrix file for parameters: " + self.params_to_string())
        matrix = self._load_matrix()
        if self.shape[0] == 0 or self.shape[1] == 0:
            logging.info("Empty matrix for parameters: " + self.params_to_string())
        (row, column, data) = matrix
        matrix = sparse.csr_matrix((data, (row, column)), shape=self.shape, dtype=float)
        return matrix

    def ref_file_available(self):
        return (self.file_name_ref is not None) and os.path.isfile(self.file_name_ref)

    def get_ref_matrix(self):
        if not self.valid:
            log.warning("No reference matrix for invalid parameters: " + self.params_to_string())
            return None
        if self.file_name_ref is None:
            raise SH.RefError("Path to reference file not specified for parameters: " + self.params_to_string())
        if not os.path.isfile(self.file_name_ref):
            raise SH.RefError("Reference basis file not found for parameters: " + self.params_to_string())
        (row, column, data) = self._load_ref_matrix()
        if len(row):
            shape=(row.pop(-1)+1,column.pop(-1)+1)
            data.pop(-1)
        else:
            shape=(0,0)
        matrix = sparse.csr_matrix((data, (row, column)), shape=shape, dtype=float)
        return matrix

    def delete_file(self):
        if os.path.isfile(self.file_name):
            os.remove(self.file_name)
