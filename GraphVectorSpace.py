from abc import ABCMeta, abstractmethod
import os
import pickle
from sage.all import  *

class GraphVectorSpace():
    __metaclass__ = ABCMeta

    def __init__(self, valid, file_name, color_counts=None, basis_on_fly=False):
        self.valid = valid
        self.file_name = file_name
        self.color_counts = color_counts
        self.basis_on_fly = basis_on_fly
        self.basis_g6 = None

        if basis_on_fly:
            if not self.basis_built():
                self.basis_g6 = self.create_basis_g6()
            else:
                self.basis_g6 = self._load_basis_g6()

    @abstractmethod
    def _generating_graphs(self):
        pass

    @abstractmethod
    def perm_sign(self, graph, perm):
        pass

    @abstractmethod
    def get_work_estimate(self):
        pass

    def create_basis_g6(self):
        outPath = self.file_name
        outDir = os.path.dirname(outPath)
        if not os.path.exists(outDir):
            os.makedirs(outDir)

        gList = self._generating_graphs()
        basisList = set()
        for G in gList:
            canonG = G.canonical_label()
            automList = G.automorphism_group().gens()
            if len(automList):
                canon6=canonG.graph6_string()
                if not canon6 in basisList:
                    if not self._has_odd_automorphisms(G, automList):
                        basisList.add(canon6)
        with open(outPath, 'wb') as f:
            pickle.dump(basisList, f)
        return basisList

    def canonical_g6(self, graph):
        canonG, permDict = graph.canonical_label(certificate=True)
        sgn = self.perm_sign(graph, [v+1 for k, v in permDict.items()])
        return (canonG.graph6_string(),sgn)

    def _has_odd_automorphisms(self, G, automList):
        for g in automList:
            if self.perm_sign(G, g.tuple()) == -1:
               return True
        return False

    def basis_built(self):
        if os.path.isfile(self.file_name):
            return True
        return False

    def _load_basis_g6(self):
        with open(self.file_name, 'rb') as f:
            basis_g6 = pickle.load(f)
        return basis_g6

    def get_basis(self, g6=False):
        if not self.valid:
            return []
        if self.basis_on_fly and self.basis_built():
            basis_g6 = self.basis_g6
        else:
            if not self.basis_built():
                raise NotBuiltError("Basis ist not built yet")
            basis_g6 = self._load_basis_g6()
        if g6:
            return basis_g6
        else:
            basisList = []
            for G in basis_g6:
                basisList.append(Graph(G))
            return basisList

    def get_dimension(self):
        if not self.valid:
            return 0
        return len(self.get_basis(g6=True))

    def delete_file(self):
        if os.path.isfile(self.file_name):
            os.remove(self.file_name)

class NotBuiltError(RuntimeError):
    pass