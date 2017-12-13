from abc import ABCMeta, abstractmethod
import os
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
    def operate_on(self,graph):
        """For G a graph returns a list of pairs (GG, x),
           such that (operator)(G) = sum x GG.
        """
        pass

    @abstractmethod
    def get_work_estimate(self):
        """Provides a rough estimate of the amount of work needed to create the operator file.
          (In arbitrary units)"""
        pass

    def create_operator_matrix(self):
        """
        Creates the matrix file that holds the operator.
        The corresponding list files for source and target
        must exist when calling this function.
        """
        try:
            domainBasis = self.domain.load_basis()
        except GVS.NotBuiltError:
            raise GVS.NotBuiltError("Cannot build operator matrix: First build basis of the domain")
        try:
            targetBasis6 = self.target.load_basis(g6=True)
        except GVS.NotBuiltError:
            raise GVS.NotBuiltError("Cannot build operator matrix: First build basis of the target")

        domainDim = len(domainBasis)
        targetDim = len(targetBasis6)

        if domainDim == 0 and targetDim == 0:
            # create empty file and return
            open(self.fileName,"w").close()
            return

        # lookup g6 -> index in target vector space
        lookup = {s: j for (j,s) in enumerate(targetBasis6)}
        matrix = []
        for (domainIndex,G) in enumerate(domainBasis):
            imageList = self._operate_on(G)
            for (GG, prefactor) in imageList:
                # canonize and look up
                GGcanon6, sgn = self.domain.canonical_g6(GG)
                imageIndex = lookup.get(GGcanon6)
                matrix.append("%d %d %d" % (domainIndex, imageIndex, sgn * prefactor))

        f = open(self.fileName, "w")
        for line in matrix:
            f.write(line + '\n')
        f.close()

    def matrix_built(self):
        if os.path.isfile(self.file_name):
            return True
        return False

    def load_operator_matrix(self):
        if not self.valid:
            return []
        if not self.matrix_built():
            raise GVS.NotBuiltError("Cannot load matrix: No matrix file")

        f = open(self.file_name, 'r')
        matrixList = f.read().splitlines()
        f.close()
        if len(matrixList)==0:
            return []
        else:
            row = []
            column = []
            data=[]
            for line in matrixList:
                (i, j, v) = map(int, line.split(" "))
                row.append(i)
                column.append(j)
                data.append(v)
            return sparse.csr_matrix(data,(row,column))

