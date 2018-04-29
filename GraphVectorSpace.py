from abc import ABCMeta, abstractmethod
from sage.all import *
import operator
import itertools
import pandas
from tqdm import tqdm
import StoreLoad as SL
import PlotCohomology
import ParallelProgress as PP
import Parameters
import Log

logger = Log.logger.getChild('graph_vector_space')


class VectorSpaceProperties(object):
    def __init__(self):
        self.dimension = None

    @classmethod
    def names(cls):
        return ['dimension']

    @staticmethod
    def sort_variables():
        return VectorSpaceProperties.names()

    def list(self):
        return [self.dimension]


class GraphVectorSpaceProperties(VectorSpaceProperties):
    def __init__(self):
        self.valid = None
        self.dimension = None

    @classmethod
    def names(cls):
        return ['valid', 'dimension']

    def list(self):
        return [self.valid, self.dimension]


class VectorSpace(object):
    __metaclass__ = ABCMeta

    def __init__(self):
        self.properties = VectorSpaceProperties()

    @abstractmethod
    def __eq__(self, other):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def get_dimension(self):
        pass

    @abstractmethod
    def get_ordered_param_dict(self):
        pass

    @abstractmethod
    def get_work_estimate(self):
        pass

    @abstractmethod
    def build_basis(self):
        pass

    @abstractmethod
    def update_properties(self):
        pass

    def get_properties(self):
        return self.properties

    def get_sort_dim(self):
        try:
            sort_dim = self.get_dimension()
        except SL.FileNotFoundError:
            sort_dim = Parameters.max_sort_value
        return sort_dim


class GraphVectorSpace(VectorSpace):
    __metaclass__ = ABCMeta

    def __init__(self):
        self.properties = GraphVectorSpaceProperties()

    @abstractmethod
    def get_type(self):
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
    def get_partition(self):
        pass

    @abstractmethod
    def is_valid(self):
        pass

    @abstractmethod
    def get_generating_graphs(self):
        """Produces a set of graphs whose isomorphism classes span the vector space. (Not necessarily freely!)"""
        pass

    @abstractmethod
    def perm_sign(self, G, p):
        """For G a graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
           Here vertex j becomes vertex p[j] in the new graph."""
        pass

    def __str__(self):
        return '<%s sub vector space with parameters: %s>' % (self.get_type(), str(self.get_ordered_param_dict()))

    def graph_to_canon_g6(self, graph):
        canonG, permDict = graph.canonical_label(partition=self.get_partition(), certificate=True)
        sgn = self.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sgn)

    def build_basis(self, progress_bar=False, ignore_existing_files=False, n_jobs=1):
        if not self.is_valid():
            return
        if not ignore_existing_files and self.exists_basis_file():
            return
        generatingList = self.get_generating_graphs()

        desc = 'Build basis: ' + str(self.get_ordered_param_dict())
        if not progress_bar:
            print(desc)
        basisSet = set()
        for G in tqdm(generatingList, desc=desc, disable=(not progress_bar)):
            if self.get_partition() is None:
                automList = G.automorphism_group().gens()
                canonG = G.canonical_label()
            else:
                automList = G.automorphism_group(partition=self.get_partition()).gens()
                canonG = G.canonical_label(partition=self.get_partition())
            if len(automList):
                canon6 = canonG.graph6_string()
                if not canon6 in basisSet:
                    if not self._has_odd_automorphisms(G, automList):
                        basisSet.add(canon6)

        #PP.parallel_progress_messaging(self._generate_basis_set, generatingList, basisSet, pbar_info=pbar_info, desc=desc)
        self._store_basis_g6(list(basisSet))

    #def _generate_basis_set(self, G, basis_set):
        #if self.get_partition() is None:
            #automList = G.automorphism_group().gens()
            #canonG = G.canonical_label()
        #else:
            #automList = G.automorphism_group(partition=self.get_partition()).gens()
            #canonG = G.canonical_label(partition=self.get_partition())
        #if len(automList):
            #canon6 = canonG.graph6_string()
            #if not canon6 in basis_set:
                #if not self._has_odd_automorphisms(G, automList):
                    #basis_set.add(canon6)

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
            logger.warn("Empty basis: %s is not valid" % str(self))
            return []
        basis_g6 = self._load_basis_g6()
        if g6:
            return basis_g6
        else:
            return [Graph(g6) for g6 in basis_g6]

    def delete_basis_file(self):
        if os.path.isfile(self.get_basis_file_path()):
            os.remove(self.get_basis_file_path())

    def update_properties(self):
        self.properties.valid = self.is_valid()
        if not self.properties.valid:
            self.properties.dimension = 0
        else:
            try:
                self.properties.dimension = self.get_dimension()
            except SL.FileNotFoundError:
                pass

    def plot_graph(self, G):
        g6 = G.graph6_string()
        path = os.path.join(self.get_plot_path(), g6 + '.png')
        SL.generate_path(path)
        P = G.plot(partition=self.get_partition(), vertex_labels=False)
        P.save(path)


class SumVectorSpace(VectorSpace):
    __metaclass__ = ABCMeta

    def __init__(self, vs_list):
        self.vs_list = vs_list
        self.is_graded = False
        try:
            if isinstance(self.vs_list[0], DegSlice):
                self.is_graded = True
        except IndexError:
            self.is_graded = True
        super(SumVectorSpace, self).__init__()

    def get_type(self):
        pass

    def get_ordered_param_range_dict(self):
        pass

    def get_ordered_param_dict(self):
        pass

    def __str__(self):
        return '<%s vector space with parameters: %s>' % (str(self.get_type()), str(self.get_ordered_param_range_dict()))

    def get_vs_list(self):
        return self.vs_list

    def get_flat_vs_list(self):
        try:
            flat_vs_list = [slice.get_vs_list() for slice in self.vs_list]
            return list(itertools.chain.from_iterable(flat_vs_list))
        except AttributeError:
            return self.vs_list

    def get_work_estimate(self):
        work_estimate = 0
        for vs in self.vs_list:
            work_estimate += vs.get_work_estimate()
        return work_estimate

    def get_dimension(self):
        dim = 0
        for vs in self.vs_list:
            dim += vs.get_dimension()
        return dim

    def update_properties(self):
        try:
            self.properties.dimension = self.get_dimension()
        except SL.FileNotFoundError:
            pass

    def contains(self, vector_space):
        for vs in self.vs_list:
            if vs == vector_space:
                return True
        return False

    def __eq__(self, other):
        if len(self.vs_list) != len(other.vs_list):
            return False
        eq_l = 0
        for (vs1, vs2) in itertools.product(self.vs_list, other.vs_list):
            if vs1 == vs2:
                eq_l += 1
        return eq_l == len(self.vs_list)

    def sort(self, key='work_estimate'):
        if key == 'work_estimate':
            self.vs_list.sort(key=operator.methodcaller('get_work_estimate'))
        elif key == 'dim':
            self.vs_list.sort(key=operator.methodcaller('get_sort_dim'))
        else:
            raise ValueError("Invalid sort key. Options: 'work_estimate', 'dim'")

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        print(' ')
        print('Build basis of %s' % str(self))
        if not (self.is_graded or isinstance(self, DegSlice)):
            self.plot_info()
        if not isinstance(self,DegSlice):
            self.sort()
        if n_jobs > 1:
            progress_bar = False
        if not self.is_graded:
            PP.parallel(self._build_single_basis, self.vs_list, n_jobs=n_jobs, progress_bar=progress_bar,
                        ignore_existing_files=ignore_existing_files)
        else:
            for vs in self.vs_list:
                self._build_single_basis(vs, progress_bar=progress_bar, ignore_existing_files=ignore_existing_files,
                                         n_jobs = n_jobs)
        self.plot_info()

    def _build_single_basis(self, vs, progress_bar=False, ignore_existing_files=True, n_jobs = 1):
        vs.build_basis(progress_bar=progress_bar, ignore_existing_files=ignore_existing_files, n_jobs = n_jobs)

    def plot_info(self):
        vsList = []
        for vs in self.vs_list:
            vs.update_properties()
            vsList.append(vs.get_ordered_param_dict().values() + vs.get_properties().list())
        try:
            vsColumns = self.vs_list[0].get_ordered_param_dict().keys() + self.vs_list[0].properties.names()
        except IndexError:
            vsColumns = []
        vsTable = pandas.DataFrame(data=vsList, columns=vsColumns)
        #vsTable.sort_values(by=VectorSpaceProperties.sort_variables(), inplace=True, na_position='last')
        PlotCohomology.display_pandas_df(vsTable)


class DegSlice(SumVectorSpace):
    def __init__(self, vs_list, deg):
        self.deg = deg
        super(DegSlice, self).__init__(vs_list)
        self.plot_info()

    def __str__(self):
        return '<degree slice of degree %d>' % self.deg

    @abstractmethod
    def get_ordered_param_dict(self):
        pass

    def get_deg(self):
        return self.deg

    def get_start_idx(self, vector_space):
        start_idx = 0
        for vs in self.vs_list:
            if vs == vector_space:
                return start_idx
            start_idx += vs.get_dimension()
        raise ValueError('vector_space should refer to a vector space of the degree slice')

    def is_complete(self):
        if len(self.vs_list) != self.deg + 1:
            return False
        for vs in self.vs_list:
            if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
                return False
        return True
