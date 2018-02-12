from abc import ABCMeta, abstractmethod
import logging
from tqdm import tqdm
from sage.all import *
import operator
import pandas
import StoreLoad as SL
import Display
import ParallelProgress as PP
import Parameters


class GraphVectorSpace():
    __metaclass__ = ABCMeta

    def __init__(self):
        self.partition = self.set_partition()
        self.valid = self.set_validity()
        self.basis_file_path = self.set_basis_file_path()
        self.plot_path = self.set_plot_path()

    @abstractmethod
    def get_params(self):
        pass

    @abstractmethod
    def get_params_string(self):
        pass

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
        try:
            dim = self.get_dimension()
        except SL.FileNotFoundError:
            dim = None
        return (self.valid, dim)

    def graph_to_canon_g6(self, graph):
        canonG, permDict = graph.canonical_label(partition=self.partition, certificate=True)
        sgn = self.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sgn)

    def build_basis(self, pbar_info, ignore_existing_files=False):
        if not self.valid:
            logging.info("Skip building basis: %s is not valid" % str(self))
            return
        if not ignore_existing_files and self.exists_basis_file():
            return
        logging.info('Build basis for graph vector space with ' + self.get_params_string())
        generatingList = self._generating_graphs()

        (progress_bar, message, idx, queue) = pbar_info
        if progress_bar:
            total = len(generatingList)
            miniters = int(total / Parameters.pbar_steps)
            desc = 'Build basis: ' + self.get_params_string()
            if message:
                queue.put((idx, 'start', total, desc))
                count = 0
            else:
                pbar = tqdm(total=total, desc=desc, miniters=miniters)
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
            if progress_bar:
                if message:
                    count += 1
                    if count % miniters == 0:
                        queue.put((idx, 'step', miniters, None))
                else:
                    pbar.update()
        if progress_bar:
            if message:
                queue.put((idx, 'stop', None, None))
            else:
                pbar.close()

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
            return int(header)
        except SL.FileNotFoundError:
            raise SL.FileNotFoundError("Dimension unknown for %s: No basis file" % str(self))

    def get_sort_value(self):
        try:
            dim = self.get_dimension()
        except SL.FileNotFoundError:
            dim = Parameters.MAX_DIMENSION
        return dim

    def _store_basis_g6(self, basisList):
        logging.info("Store basis in file: %s" % str(self.basis_file_path))
        basisList.insert(0, str(len(basisList)))
        SL.store_string_list(basisList, self.basis_file_path)

    def _load_basis_g6(self):
        if not self.exists_basis_file():
            raise SL.FileNotFoundError("Cannot load basis, No Basis file found for %s: " % str(self))
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

    def plot_graph(self, G):
        g6 = G.graph6_string()
        path = os.path.join(self.plot_path, g6 + '.png')
        SL.generate_path(path)
        P = G.plot(partition=self.partition, vertex_labels=False)
        P.save(path)


class VectorSpaceCollection:
    __metaclass__ = ABCMeta

    def __init__(self, vs_list):
        self.vs_list = vs_list

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_params_range(self):
        pass

    @abstractmethod
    def get_params_names(self):
        pass

    @classmethod
    def sum(cls, vs_collections_1, vs_collections_2):
        return cls(vs_collections_1.vs_list + vs_collections_2.vs_list)

    def sort(self, work_estimate=True):
        if work_estimate:
            self.vs_list.sort(key=operator.methodcaller('get_work_estimate'))
        else:
            self.vs_list.sort(key=operator.methodcaller('get_sort_value'))

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.plot_info()
        self.sort()
        PP.parallel_individual_progress(self._build_single_basis, self.vs_list, n_jobs=n_jobs,
                                        progress_bar=progress_bar, ignore_existing_files=ignore_existing_files)

    def _build_single_basis(self, vs, pbar_info, ignore_existing_files=True):
        vs.build_basis(pbar_info, ignore_existing_files=ignore_existing_files)

    def plot_info(self):
        vsList = []
        for vs in self.vs_list:
            vsList.append(list(vs.get_params()) + list(vs.get_info()))
        vsColumns = list(self.get_params_names()) + ['valid', 'dimension']
        vsTable = pandas.DataFrame(data=vsList, columns=vsColumns)
        vsTable.sort_values(by=['valid', 'dimension'], inplace=True, na_position='last')
        Display.display_pandas_df(vsTable)
