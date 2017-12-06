from abc import ABCMeta, abstractmethod
import os
from sage.all import  *

class GraphVectorSpace():
    __metaclass__ = ABCMeta
    @abstractmethod
    def get_file_name(self):
        pass

    @abstractmethod
    def get_svg_dir(self):
        pass

    @abstractmethod
    def is_valid(self):
        pass

    @abstractmethod
    def get_generating_graphs(self):
        pass

    @abstractmethod
    def get_color_counts(self):
        pass

    @abstractmethod
    def get_perm_sign(self, graph, perm):
        pass

    @abstractmethod
    def get_work_estimate(self):
        pass

    def hasOddAutomorphisms(self, G, automList):
        for p in automList:
            if self.get_perm_sign(G, p) == -1:
                return True
        return False

    def createListFile(self):
        outPath = self.get_file_name()
        outDir = os.path.dirname(outPath)
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        colorData = self.get_color_counts()
        gL = self.get_generating_graphs()

        gS = set()
        for G in gL:
            canonG, automList = self.get_canon_and_automorphisms(G, colorData)
            canon6=canonG.graph6_string()
            if not canon6 in gS:
                if not self.hasOddAutomorphisms(G, automList):
                    gS.add(canon6)

        f = open(outPath, 'w')
        for g6 in gS:
            f.write(g6 + '\n')
        f.close()

    def getDimension(self):
        if not self.is_valid():
            return 0
        fileName = self.get_file_name()
        if not os.path.isfile(fileName):
            return -1
        f = open(fileName, 'r')
        dimension=0
        for line in f:
            dimension += 1
        f.close()
        return dimension