from abc import ABCMeta, abstractmethod
import os
from sage.all import  *
import GraphVectorSpace as GVS

class GraphOperator():
    __metaclass__ = ABCMeta
    @abstractmethod
    def file_name(self):
        """Retrieve the file name (and path) of the file storing the matrix."""
        pass

    @abstractmethod
    def domain(self):
        """Returns the GraphVectorSpace on which the operator acts."""
        pass

    @abstractmethod
    def target(self):
        """Returns the GraphVectorSpace in which the operator takes values."""
        pass

    @abstractmethod
    def _operate_on(self,graph):
        """For G a graph returns a list of pairs (GG, x),
           such that (operator)(G) = sum x GG.
        """
        pass

    @abstractmethod
    def work_estimate(self):
        """Provides a rough estimate of the amount of work needed to create the operator file.
          (In arbitrary units)"""
        pass

    def valid(self):
        return self.domain().valid() and self.target().valid()

    def create_operator_matrix(self):
        """
        Creates the matrix file that holds the operator.
        The corresponding list files for source and target
        must exist when calling this function.
        """
        fileName = self.file_name()
        domain = self.domain()
        target = self.target()

        try:
            domainBasis = domain.basis()
        except GVS.NotBuiltError:
            raise GVS.NotBuiltError("Cannot build operator matrix: First build basis of the domain")
        try:
            targetBasis6 = target.basis(g6=True)
        except GVS.NotBuiltError:
            raise GVS.NotBuiltError("Cannot build operator matrix: First build basis of the target")

        domainDim = len(domainBasis)
        targetDim = len(targetBasis6)

        if domainDim == 0 and targetDim == 0:
            # create empty file and return
            open(fileName,"w").close()
            return

        color_counts = domain.color_counts()
        matrix = []

        # lookup g6 -> index in target vector space
        lookup = {s: j for (j,s) in enumerate(targetBasis6)}

        f = open(fileName,"w")
        for (domainIndex,G) in enumerate(domainBasis):
            imageList = self.operate_on(G)
            for (GG, prefactor) in imageList:
                # canonize and look up
                GGcanon6, sgn = domain.canonical_g6(GG)
                imageIndex = lookup.get(GGcanon6)
                if imageIndex:
                    f.write("%d %d %d\n" % (domainIndex, imageIndex, sgn * prefactor))
          # write matrix size