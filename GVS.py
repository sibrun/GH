from abc import ABC, abstractmethod
import os

class GraphVectorSpace(ABC):
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
    def get_dot(self, graph):
        pass

    @abstractmethod
    def get_work_estimate(self):
        pass

    def hasOddAutomorphisms(self, graph, automList):
        for p in automList
            if self.get_perm_sign(graph, p) == -1:
                return True
        return False

    def createListFile(self):

        outPath = self.get_file_name()
        outDir = os.path.dirname(outPath)
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        colorData = self.get_color_counts()
        graphList = self.get_generating_graphs()

        graphSet = set()
        for graph in graphList:
            canonGraph, automList = self.get_canon_and_automorphisms(graph, colorData)
            canon6 = self.to_string(canonGraph)
            if not canon6 in graphSet:
                if not self.hasOddAutomorphisms(graph, automList):
                    graphSet.add(graph)

        f = open(outPath, 'w')
        for g6 in graphSet:
            f.write(g6 + '\n')
        f.close()

    def getDimension(self):
        if not self.is_valid():
            return 0
        fileName = self.get_file_name()
        if not os.path.isfile(fileName):
            return -1
        f = open(outPath, 'r')
        dimension=0
        for line in f:
            dimension += 1
        return dimension