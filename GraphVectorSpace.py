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
        self.valid = self._set_validity()
        self.file_path = self._set_file_path()
        self.work_estimate = self._set_work_estimate()
        self.file_path_ref = self._set_file_path(ref=True)
        self.header_ref = header_ref
        self.dimension = None
        self.color_counts = color_counts

    @abstractmethod
    def _set_file_path(self, ref=False):
        pass

    @abstractmethod
    def _set_validity(self,):
        pass

    @abstractmethod
    def _set_work_estimate(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    @abstractmethod
    def _generating_graphs(self):
        pass

    @abstractmethod
    def perm_sign(self, graph, perm):
        pass

    def get_info(self):
        validity = "valid" if self.valid else "not valid"
        built = "basis built" if self.basis_built() else "basis not built"
        dimension = "dimension unknown"
        if self.basis_built():
            dimension = "dimension = none" if self.dimension is None else "dim = %d" % self.dimension
        return "%s, %s, %s" % (validity, built, dimension)

    def graph_to_canon_g6(self, graph):
        canonG, permDict = graph.canonical_label(certificate=True)
        sgn = self.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sgn)

    def g6_to_canon_g6(self, graph6):
        graph = Graph(graph6)
        return graph.canonical_label().graph6_string()

    def build_basis(self):
        if not self.valid:
            logging.info("Skip building basis: %s is not valid" % str(self))
            return
        if self.basis_built():
            self.load_info_from_file()
            return
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
        self._store_basis_g6(list(basisSet))
        logging.info("Basis built for %s" % str(self))

    def _has_odd_automorphisms(self, G, automList):
        for g in automList:
            if self.perm_sign(G, g.tuple()) == -1:
               return True
        return False

    def basis_built(self):
        if os.path.isfile(self.file_path):
            return True
        return False

    def load_info_from_file(self):
        try:
            header = SH.load_header(self.file_path)
        except SH.FileNotExistingError:
            logging.warn("Cannot load infos from file %s" % str(self.file_path))
            return
        self.dimension = int(header)
        logging.info("Infos loaded from file %s" % str(self.file_path))

    def _store_basis_g6(self, basis_g6):
        logging.info("Store basis in file: %s" % str(self.file_path))
        SH.store_string_list(basis_g6, self.file_path, header=str(self.dimension))

    def _load_basis_g6(self):
        logging.info("Load basis from file: %s" % str(self.file_path))
        (header, basis_g6) = SH.load_string_list(self.file_path, header=True)
        self.dimension = int(header)
        if len(basis_g6) != self.dimension:
            raise ValueError("Basis read from file %s has wrong dimension" % str(self.file_path))
        return basis_g6

    def _load_ref_basis_g6(self):
        logging.info("Load basis from reference file: %s" % str(self.file_path_ref))
        data = SH.load_string_list(self.file_path_ref, header=self.header_ref)
        if self.header_ref:
            (header, basis_g6) = data
            if len(basis_g6) != int(header):
                raise ValueError("Basis read from file %s has wrong dimension" % str(self.file_path_ref))
            return basis_g6
        else:
            return data

    def get_basis(self, g6=True):
        if not self.valid:
            logging.warn("Empty basis: %s is not valid" % str(self))
            return []
        if not self.basis_built():
            raise SH.NotBuiltError("Cannot load basis, No Basis file found for %s: " % str(self))
        basis_g6 = self._load_basis_g6()
        logging.info("Get basis of %s with dimension %d" % (str(self), self.dimension))
        if g6:
            return basis_g6
        else:
            basis = []
            for G in basis_g6:
                basis.append(Graph(G))
            return basis

    def ref_file_available(self):
        return (self.file_path_ref is not None) and os.path.isfile(self.file_path_ref)

    def get_ref_basis_g6(self):
        if self.file_path_ref is None:
            raise SH.RefError("%s: Path to reference file not specified" % str(self))
        if not os.path.isfile(self.file_path_ref):
            raise SH.RefError("%s: Reference basis file not found" % str(self))
        basis_g6 = self._load_ref_basis_g6()
        basis_g6_canon = []
        for G6 in basis_g6:
            canonG6 = self.g6_to_canon_g6(G6)
            basis_g6_canon.append(canonG6)
        return basis_g6_canon

    def delete_file(self):
        if os.path.isfile(self.file_path):
            os.remove(self.file_path)
