from abc import ABCMeta, abstractmethod
import logging
from sage.all import *
import StoreLoad as SL


class GraphVectorSpace():
    __metaclass__ = ABCMeta

    def __init__(self):
        self.partition = self.set_partition()
        self.valid = self.set_validity()
        self.basis_file_path = self.set_basis_file_path()
        self.plot_path = self.set_plot_path()

    @abstractmethod
    def set_partition(self):
        pass

    @abstractmethod
    def set_basis_file_path(self):
        pass

    @abstractmethod
    def set_plot_path(self):
        pass

    @abstractmethod
    def get_ref_basis_file_path(self):
        pass

    @abstractmethod
    def set_validity(self, ):
        pass

    @abstractmethod
    def get_work_estimate(self):
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
    def perm_sign(self, G, p):
        pass

    def get_info(self):
        if not self.valid:
            return "not valid"
        dim = "dimension unknown"
        if self.exists_basis_file():
            dim = "dimension: %d" % self.get_dimension()
        return "valid, %s" % dim

    def graph_to_canon_g6(self, graph):
        canonG, permDict = graph.canonical_label(certificate=True)
        sign = self.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sign)

    def build_basis(self, ignore_existing_file=False):
        if not self.valid:
            logging.info("Skip building basis: %s is not valid" % str(self))
            return
        if not ignore_existing_file and self.exists_basis_file():
            return
        generatingList = self._generating_graphs()
        basisSet = set()
        for G in generatingList:
            if self.partition is None:
                automList = G.automorphism_group().gens()
                canonG = G.canonical_label()
            else:
                automList = G.automorphism_group(partition=self.partition).gens()
                canonG = G.canonical_label(partition=self.partition)
            if len(automList):
                canon6=canonG.graph6_string()
                if not canon6 in basisSet:
                    if not self._has_odd_automorphisms(G, automList):
                        basisSet.add(canon6)
        self._store_basis_g6(list(basisSet))
        logging.info("Basis built for %s" % str(self))

    def _has_odd_automorphisms(self, G, automList):
        for g in automList:
            if self.perm_sign(G, g.tuple()) == -1:
               return True
        return False

    def exists_basis_file(self):
        return os.path.isfile(self.basis_file_path)

    def get_dimension(self):
        if not self.valid:
            return 0
        try:
            header = SL.load_line(self.basis_file_path)
        except SL.FileNotExistingError:
            raise SL.NotBuiltError("Cannot load header from file %s: Build basis first" % str(self.basis_file_path))
        return int(header)

    def _store_basis_g6(self, basisList):
        logging.info("Store basis in file: %s" % str(self.basis_file_path))
        basisList.insert(0, str(len(basisList)))
        SL.store_string_list(basisList, self.basis_file_path)

    def _load_basis_g6(self):
        if not self.exists_basis_file():
            raise SL.NotBuiltError("Cannot load basis, No Basis file found for %s: " % str(self))
        logging.info("Load basis from file: %s" % str(self.basis_file_path))
        basisList = SL.load_string_list(self.basis_file_path)
        dim = int(basisList.pop(0))
        if len(basisList) != dim:
            raise ValueError("Basis read from file %s has wrong dimension" % str(self.basis_file_path))
        return basisList

    def get_basis(self, g6=True):
        if not self.valid:
            logging.warn("Empty basis: %s is not valid" % str(self))
            return []
        basis_g6 = self._load_basis_g6()
        logging.info("Get basis of %s with dimension %d" % (str(self), len(basis_g6)))
        if g6:
            return basis_g6
        else:
            return [Graph(g6) for g6 in basis_g6]

    def delete_basis_file(self):
        if os.path.isfile(self.basis_file_path):
            os.remove(self.basis_file_path)
