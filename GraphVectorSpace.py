from abc import ABCMeta, abstractmethod
import os
from sage.all import  *

class GraphVectorSpace():
    __metaclass__ = ABCMeta
    @abstractmethod
    def file_name(self):
        pass

    @abstractmethod
    def svg_dir(self):
        pass

    @abstractmethod
    def valid(self):
        pass

    @abstractmethod
    def _generating_graphs(self):
        pass

    @abstractmethod
    def color_counts(self):
        pass

    @abstractmethod
    def _perm_sign(self, graph, perm):
        pass

    @abstractmethod
    def canonical(self, graph, colorData):
        pass

    @abstractmethod
    def work_estimate(self):
        pass

    def create_basis(self):
        outPath = self.file_name()
        outDir = os.path.dirname(outPath)
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        gL = self._generating_graphs()
        gS = set()
        for G in gL:
            canonG = G.canonical_label()
            automList = G.automorphism_group().gens()
            if len(automList):
                canon6=canonG.graph6_string()
                if not canon6 in gS:
                    if not self._has_odd_automorphisms(G, automList):
                        gS.add(canon6)

        f = open(outPath, 'w')
        for g6 in gS:
            f.write(g6 + '\n')
        f.close()

    def _has_odd_automorphisms(self, G, automList):
        for g in automList:
            if self._perm_sign(G, g.tuple()) == -1:
               return True
        return False

    def basis_built(self):
        if os.path.isfile(self.file_ame()):
            return True
        return False

    def dimension(self):
        if not self.valid():
            return 0
        fileName=self.file_name()
        if not os.path.isfile(fileName):
            raise NotBuiltError("Bsais ist not built yet")
        f = open(fileName, 'r')
        dimension=0
        for line in f:
            dimension += 1
        f.close()
        return dimension

    def basis(self):
        if not self.valid():
            return []
        fileName = self.file_name()
        if not os.path.isfile(fileName):
            raise NotBuiltError("Bsais ist not built yet")
        f = open(fileName, 'r')
        basisList=[]
        for line in f:
            basisList.append(Graph(line))
        f.close()
        return basisList

class NotBuiltError(RuntimeError):
    pass