from abc import ABCMeta, abstractmethod
import os
import pickle
import scipy.sparse as sparse
from sage.all import  *
import GraphVectorSpace as GVS

class GraphOperator():
    __metaclass__ = ABCMeta

    def __init__(self, file_name, domain, target):
        self.file_name = file_name
        self.domain = domain
        self.target = target
        self.valid = self.domain.valid and self.target.valid

    @abstractmethod
    def _operate_on(self,graph):
        """For G a graph returns a list of pairs (GG, x),
           such that (operator)(G) = sum x GG.
        """
        pass

    @abstractmethod
    def get_work_estimate(self):
        """Provides a rough estimate of the amount of work needed to create the operator file.
          (In arbitrary units)"""
        pass

    def create_domain_basis(self):
        self.domain.create_basis()

    def create_target_basis(self):
        self.target.create_basis()

    def create_operator_matrix(self):
        """
        Creates the matrix file that holds the operator.
        The corresponding list files for source and target
        must exist when calling this function.
        """
        if not self.matrix_built():
            try:
                domainBasis = self.domain.get_basis(g6=False)
            except GVS.NotBuiltError:
                raise GVS.NotBuiltError("Cannot build operator matrix: First build basis of the domain")
            try:
                targetBasis6 = self.target.get_basis(g6=True)
            except GVS.NotBuiltError:
                raise GVS.NotBuiltError("Cannot build operator matrix: First build basis of the target")

            domainDim = len(domainBasis)
            targetDim = len(targetBasis6)

            if domainDim == 0 or targetDim == 0:
                # create empty file and return
                open(self.fileName,"w").close()
                return
            # lookup g6 -> index in target vector space
            lookup = {s: j for (j,s) in enumerate(targetBasis6)}
            print(lookup)
            matrix = []
            for (domainIndex,G) in enumerate(domainBasis):
                imageList = self._operate_on(G)
                for (GG, prefactor) in imageList:
                    # canonize and look up
                    GGcanon6, sgn = self.domain.canonical_g6(GG)
                    #print(GGcanon6)
                    imageIndex = lookup.get(GGcanon6)
                    matrix.append((domainIndex, imageIndex, sgn * prefactor))

            #print(matrix)
            self._store_operator_matrix(matrix)

    def matrix_built(self):
        if os.path.isfile(self.file_name):
            return True
        return False

    def _store_operator_matrix(self, matrixList):
        with open(self.file_name, 'w') as f:
            for line in matrixList:
                (i,j,v) = line
                line_string = "%d %d %d" % (i,j,v)
                f.write(line_string + '\n')

    def _load_operator_matrix(self):
        with open(self.file_name, 'r') as f:
            matrixList_string = f.read().splitlines()
        matrixList=[]
        for line in matrixList_string:
            (i, j, v) = map(int, line.split(" "))
            matrixList.append((i, j, v))
        return matrixList

    def get_operator_matrix(self):
        if not self.valid:
            return []
        if not self.matrix_built():
            raise GVS.NotBuiltError("Cannot load matrix: No matrix file")

        matrixList = self._load_operator_matrix()
        if len(matrixList)==0:
            return []

        row = []
        column = []
        data=[]
        for (i, j, v) in matrixList:
            row.append(i)
            column.append(j)
            data.append(v)
        return sparse.csr_matrix(data,(row,column))

    def delete_file(self):
        if os.path.isfile(self.file_name):
            os.remove(self.file_name)
