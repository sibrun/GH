from abc import ABCMeta, abstractmethod
import os
from sage.all import  *

class GraphVectorSpace():
    __metaclass__ = ABCMeta

    def __init__(self, valid, file_name, color_counts=None):
        self.valid = valid
        self.file_name = file_name
        self.color_counts = color_counts

    @abstractmethod
    def _generating_graphs(self):
        pass

    @abstractmethod
    def perm_sign(self, graph, perm):
        pass

    @abstractmethod
    def get_work_estimate(self):
        pass

    def create_basis(self):
        if not self.basis_built():
            generatingList = self._generating_graphs()
            basisSet = set()
            for G in generatingList:
                canonG = G.canonical_label()
                automList = G.automorphism_group().gens()
                if len(automList):
                    canon6=canonG.graph6_string()
                    if not canon6 in basisSet:
                        if not self._has_odd_automorphisms(G, automList):
                            basisSet.add(canon6)
            self._store_basis_g6(list(basisSet))

    def canonical_g6(self, graph):
        canonG, permDict = graph.canonical_label(certificate=True)
        sgn = self.perm_sign(graph, permDict.values())
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

    def _store_basis_g6(self, basis_g6):
        outDir = os.path.dirname(self.file_name)
        if not os.path.exists(outDir):
            os.makedirs(outDir)
        with open(self.file_name, 'w') as f:
            for g6 in basis_g6:
                f.write(g6 + '\n')

    def _load_basis_g6(self):
        with open(self.file_name, 'r') as f:
            basis_g6 = f.read().splitlines()
        return basis_g6

    def get_basis(self, g6=True):
        if not self.valid:
            return []
        if not self.basis_built():
            raise NotBuiltError("Basis ist not built yet")
        basis_g6 = self._load_basis_g6()
        if g6:
            return basis_g6
        else:
            basis = []
            for G in basis_g6:
                basis.append(Graph(G))
            return basis

    def get_dimension(self):
        if not self.valid:
            return 0
        return len(self.get_basis(g6=True))

    def delete_file(self):
        if os.path.isfile(self.file_name):
            os.remove(self.file_name)

class NotBuiltError(RuntimeError):
    pass