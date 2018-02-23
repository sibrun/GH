from abc import ABCMeta, abstractmethod
import logging
from sage.all import *
import operator
import itertools
import pandas
import StoreLoad as SL
import Display
import ParallelProgress as PP
import Parameters


class GraphVectorSpace(object):
    __metaclass__ = ABCMeta

    @abstractmethod
    def get_params_dict(self):
        pass

    @abstractmethod
    def get_partition(self):
        pass

    @abstractmethod
    def get_basis_file_path(self):
        pass

    @abstractmethod
    def get_plot_path(self):
        pass

    @abstractmethod
    def get_ref_basis_file_path(self):
        pass

    @abstractmethod
    def is_valid(self, ):
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
    def generating_graphs(self):
        pass

    @abstractmethod
    def perm_sign(self, G, p):
        pass

    def get_info_dict(self):
        try:
            dim = self.get_dimension()
        except SL.FileNotFoundError:
            dim = None
        return {'valid': self.is_valid(), 'dimension': dim}

    def graph_to_canon_g6(self, graph):
        canonG, permDict = graph.canonical_label(partition=self.get_partition(), certificate=True)
        sgn = self.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sgn)

    def build_basis(self, pbar_info=False, ignore_existing_files=False):
        if not self.is_valid():
            return
        if not ignore_existing_files and self.exists_basis_file():
            return
        generatingList = self.generating_graphs()

        desc = 'Build basis: ' + str(self.get_params_dict())
        basisSet = set()
        PP.parallel_progress_messaging(self._generate_basis_set, generatingList, basisSet,
                                       pbar_info=pbar_info, desc=desc)
        self._store_basis_g6(list(basisSet))

    def _generate_basis_set(self, G, basis_set):
        if self.get_partition() is None:
            automList = G.automorphism_group().gens()
            canonG = G.canonical_label()
        else:
            automList = G.automorphism_group(partition=self.get_partition()).gens()
            canonG = G.canonical_label(partition=self.get_partition())
        if len(automList):
            canon6 = canonG.graph6_string()
            if not canon6 in basis_set:
                if not self._has_odd_automorphisms(G, automList):
                    basis_set.add(canon6)

    def _has_odd_automorphisms(self, G, automList):
        for g in automList:
            if self.perm_sign(G, g.tuple()) == -1:
               return True
        return False

    def exists_basis_file(self):
        return os.path.isfile(self.get_basis_file_path())

    def get_dimension(self):
        if not self.is_valid():
            return 0
        try:
            header = SL.load_line(self.get_basis_file_path())
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
        basisList.insert(0, str(len(basisList)))
        SL.store_string_list(basisList, self.get_basis_file_path())

    def _load_basis_g6(self):
        if not self.exists_basis_file():
            raise SL.FileNotFoundError("Cannot load basis, No basis file found for %s: " % str(self))
        basisList = SL.load_string_list(self.get_basis_file_path())
        dim = int(basisList.pop(0))
        if len(basisList) != dim:
            raise ValueError("Basis read from file %s has wrong dimension" % str(self.get_basis_file_path()))
        return basisList

    def get_basis(self, g6=True):
        if not self.is_valid():
            logging.warn("Empty basis: %s is not valid" % str(self))
            return []
        basis_g6 = self._load_basis_g6()
        if g6:
            return basis_g6
        else:
            return [Graph(g6) for g6 in basis_g6]

    def delete_basis_file(self):
        if os.path.isfile(self.get_basis_file_path()):
            os.remove(self.get_basis_file_path())

    def plot_graph(self, G):
        g6 = G.graph6_string()
        path = os.path.join(self.get_plot_path(), g6 + '.png')
        SL.generate_path(path)
        P = G.plot(partition=self.get_partition(), vertex_labels=False)
        P.save(path)


class VectorSpace(object):
    __metaclass__ = ABCMeta

    def __init__(self, vs_list):
        self.vs_list = vs_list

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_params_range_dict(self):
        pass

    def get_vs_list(self):
        return self.vs_list

    def __eq__(self, other):
        if len(self.vs_list) != len(other.vs_list):
            return False
        eq_l = 0
        for (vs1, vs2) in itertools.product(self.vs_list, other.vs_list):
            if vs1 == vs2:
                eq_l += 1
        return eq_l == len(self.vs_list)

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

    def _build_single_basis(self, vs, pbar_info=False, ignore_existing_files=True):
        vs.build_basis(pbar_info=pbar_info, ignore_existing_files=ignore_existing_files)

    def plot_info(self):
        vsList = []
        for vs in self.vs_list:
            vsList.append(vs.get_params_dict().values() + vs.get_info_dict().values())
        vsColumns = self.get_params_range_dict().keys() + ['valid', 'dimension']
        vsTable = pandas.DataFrame(data=vsList, columns=vsColumns)
        vsTable.sort_values(by=['valid', 'dimension'], inplace=True, na_position='last')
        Display.display_pandas_df(vsTable)


class URBiGradedVectorSpace(VectorSpace):
    __metaclass__ = ABCMeta

    def __init__(self, vs_list):
        super(URBiGradedVectorSpace, self).__init__(vs_list)
        self.graded_vs_list = []
        self.build_grading()

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_params_range_dict(self):
        pass

    @abstractmethod
    def get_deg_params(self):
        pass

    def build_grading_dict(self):
        deg_params = self.get_deg_params()
        grading_dict = {}
        if len(deg_params) != 2:
            raise ValueError('Bigraded vector space needs two degree parameters')
        for vs in self.vs_list:
            try:
                tot_deg = getattr(vs, deg_params[0]) + getattr(vs, deg_params[1])
            except AttributeError:
                raise AttributeError('Vector space %s does not have the expected degree parameters %s, %s'
                                     % (str(vs), deg_params[0], deg_params[1]))
            deg_slice = grading_dict.get(tot_deg)
            if deg_slice is None:
                grading_dict.update({tot_deg, DegSliceURBiGradedVS(tot_deg, deg_params[0])})
            else:
                grading_dict[tot_deg].append(vs)
        graded_vs_list = grading_dict.values()
        for degSlice in graded_vs_list:
            degSlice.test_completeness()
        self.graded_vs_list = graded_vs_list


class DegSliceURBiGradedVS(object):
    def __init__(self, tot_deg, idx_deg):
        self.vs_list = [None] * tot_deg + 1
        self.idx_deg = idx_deg

    def get_dimension(self):
        dim = 0
        for vs in self.vs_list:
            dim += vs.get_dimension()

    def append(self, vs):
        try:
            idx = vs.getattr(self.idx_deg)
        except AttributeError:
            raise AttributeError('Vector space %s does not hav the expectet degree %s' % (str(vs), self.idx_deg))
        try:
            self.vs_list[idx] = vs
        except IndexError:
            raise ValueError('Vector space %s does not have a valid index in the degree slice' % str(vs))

    def test_completeness(self):
        for vs in self.vs_list:
            if vs is None or (vs.is_valid() and vs.exists_bsis_file()):
                raise AssertionError('Degree slice %s is not complete: Build complete basis first' % str(self))
        return True
