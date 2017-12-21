from abc import ABCMeta, abstractmethod
import os
import logging
import pickle
from sage.all import *
import Shared as SH

reload(SH)

class GraphVectorSpace():
    __metaclass__ = ABCMeta

    def __init__(self, header_ref=False, color_counts=None):
        self.valid = self._set_validity
        self.file_name = self._set_file_name()
        self.work_estimate = self._set_work_estimate()
        self.file_name_ref = self._set_file_name(ref=True)
        self.header_ref = header_ref
        self.dimension = None
        self.color_counts = color_counts

    @abstractmethod
    def _set_file_name(self, ref=False):
        pass

    @abstractmethod
    def _set_validity(self,):
        pass

    @abstractmethod
    def _set_work_estimate(self):
        pass

    @abstractmethod
    def _generating_graphs(self):
        pass

    @abstractmethod
    def perm_sign(self, graph, perm):
        pass

    @abstractmethod
    def params_to_string(self):
        pass

    def graph_to_canon_g6(self, graph):
        canonG, permDict = graph.canonical_label(certificate=True)
        sgn = self.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sgn)

    def g6_to_canon_g6(self, graph6):
        graph = Graph(graph6)
        return graph.canonical_label().graph6_string()

    def build_basis(self):
        if not self.valid:
            logging.info("Skip building vector space basis for invalid parameters: " + self.params_to_string())
            return
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
            self.dimension = len(basisSet)
            logging.info("Vector space basis built for parameters: " + self.params_to_string())
            self._store_basis_g6(list(basisSet))

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
        logging.info("Store basis in file: " + self.file_name)
        out_dir = os.path.dirname(self.file_name)
        if not os.path.exists(out_dir):
            os.makedirs(out_dir)
        with open(self.file_name, 'w') as f:
            f.write(str(self.dimension) + '\n')
            if self.dimension:
                for g6 in basis_g6:
                    f.write(g6 + '\n')

    def _load_basis_g6(self):
        logging.info("Access basis file: " + self.file_name)
        with open(self.file_name, 'r') as f:
            self.dimension = int(f.readline())
            if self.dimension == 0:
                return []
            basis_g6 = f.read().splitlines()
            if len(basis_g6) != self.dimension:
                raise ValueError("Basis read from file has wrong dimension")
            return basis_g6

    def _load_ref_basis_g6(self):
        logging.info("Access reference file: " + self.file_name_ref)
        with open(self.file_name_ref, 'r') as f:
            if self.header_ref:
                dim = int(f.readline())
            return f.read().splitlines()

    def get_basis(self, g6=True):
        if not self.valid:
            logging.warn("Empty basis for invalid parameters: " + self.params_to_string())
            return []
        if not self.basis_built():
            raise SH.NotBuiltError("Cannot load basis, No Basis file for parameters: " + self.params_to_string())
        basis_g6 = self._load_basis_g6()
        if self.dimension == 0:
            logging.info("Empty basis for parameters: " + self.params_to_string())
        if g6:
            return basis_g6
        else:
            basis = []
            for G in basis_g6:
                basis.append(Graph(G))
            return basis

    def ref_file_available(self):
        return (self.file_name_ref is not None) and os.path.isfile(self.file_name_ref)

    def get_ref_basis_g6(self):
        if self.file_name_ref is None:
            raise SH.RefError("Path to reference file not specified for parameters: " + self.params_to_string())
        if not os.path.isfile(self.file_name_ref):
            raise SH.RefError("Reference basis file not found for parameters: " + self.params_to_string())
        basis_g6 = self._load_ref_basis_g6()
        basis_g6_canon = []
        for G6 in basis_g6:
            canonG6 = self.g6_to_canon_g6(G6)
            basis_g6_canon.append(canonG6)
        return basis_g6_canon

    def delete_file(self):
        if os.path.isfile(self.file_name):
            os.remove(self.file_name)
