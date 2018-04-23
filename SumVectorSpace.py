from abc import ABCMeta, abstractmethod
from sage.all import *
import operator
import itertools
import pandas
from tqdm import tqdm
import StoreLoad as SL
import Display
import ParallelProgress as PP
import Parameters
import Log

logger = Log.logger.getChild('graph_vector_space')


class VectorSpaceProperties(object):
    def __init__(self):
        self.valid = None
        self.dimension = None

    @staticmethod
    def names():
        return ['valid', 'dimension']

    @staticmethod
    def sort_variables():
        return VectorSpaceProperties.names()

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
    def is_valid(self):
        pass

    @abstractmethod
    def get_dimension(self):
        pass

    def get_properties(self):
        return self.properties


class GraphVectorSpace(VectorSpace):
    __metaclass__ = ABCMeta

    def __init__(self):
        super(GraphVectorSpace, self).__init__()

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_ordered_param_dict(self):
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
    def get_work_estimate(self):
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

    def build_basis(self, progress_bar=False, ignore_existing_files=False):
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

    def get_sort_dim(self):
        try:
            sort_dim = self.get_dimension()
        except SL.FileNotFoundError:
            sort_dim = Parameters.max_sort_value
        return sort_dim

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
        elif self.exists_basis_file():
            self.properties.dimension = self.get_dimension()

    def get_properties(self):
        return self.properties

    def plot_graph(self, G):
        g6 = G.graph6_string()
        path = os.path.join(self.get_plot_path(), g6 + '.png')
        SL.generate_path(path)
        P = G.plot(partition=self.get_partition(), vertex_labels=False)
        P.save(path)


class SumVectorSpace(object):
    __metaclass__ = ABCMeta

    def __init__(self, vs_list):
        self.vs_list = vs_list

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_ordered_param_range_dict(self):
        pass

    def __str__(self):
        return '<%s vector space with parameters: %s>' % (self.get_type(), str(self.get_ordered_param_range_dict()))

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
        self.plot_info()
        self.sort()
        if n_jobs > 1:
            progress_bar = False
        PP.parallel(self._build_single_basis, self.vs_list, n_jobs=n_jobs, progress_bar=progress_bar,
                    ignore_existing_files=ignore_existing_files)
        self.plot_info()

    def _build_single_basis(self, vs, progress_bar=False, ignore_existing_files=True):
        vs.build_basis(progress_bar=progress_bar, ignore_existing_files=ignore_existing_files)

    def plot_info(self):
        vsList = []
        for vs in self.vs_list:
            vs.update_properties()
            vsList.append(vs.get_ordered_param_dict().values() + vs.get_properties().list())
        vsColumns = self.get_ordered_param_range_dict().keys() + VectorSpaceProperties.names()
        vsTable = pandas.DataFrame(data=vsList, columns=vsColumns)
        #vsTable.sort_values(by=VectorSpaceProperties.sort_variables(), inplace=True, na_position='last')
        Display.display_pandas_df(vsTable)


class DegSlice(VectorSpace):
    def __init__(self, deg):
        self.deg = deg
        self.vs_dict = dict()
        super(DegSlice, self).__init__()

    def __str__(self):
        return '<degree slice of degree %d>' % self.deg

    def __eq__(self, other):
        return self.vs_dict.items() == other.vs_dict.items()

    def is_valid(self):
        all_not_valid = True
        for vs in self.vs_dict.keys():
            if vs.is_valid():
                all_not_valid = False
        return not all_not_valid

    def get_deg(self):
        return self.deg

    def get_dimension(self):
        dim = 0
        for vs in self.vs_dict.keys():
            dim += vs.get_dimension()

    def get_vs_list(self):
        return self.vs_dict.keys()

    def get_start_idx(self, sub_vector_space):
        if self.vs_dict.get(sub_vector_space) is None:
            raise ValueError('sub_vector_space needs to refer on a vector space of the degree slice')
        vs_list = self.get_vs_list()
        start_idx = 0
        for vs in vs_list:
            if vs == sub_vector_space:
                return start_idx
            start_idx += vs.get_dimension()

    def append(self, vs, idx):
        self.vs_dict.update({vs: idx})

    def is_complete(self):
        for idx in range(0, self.deg + 2):
            vs = self.vs_dict.get(idx)
            if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
                return False
        return True

class BiGrading(object):
    def __init__(self, vector_space):
        self.vector_space = vector_space
        self.grading_dict = dict()
        self.build_grading()

    @abstractmethod
    def get_degs(self, graph_vs):
        pass

    def get_deg_slices(self):
        return self.grading_dict.values()

    def build_grading(self):
        for vs in self.vector_space.get_vs_list():
            (deg1, deg2) = self.get_degs(vs)
            deg_tot = deg1 + deg2
            idx = deg1
            deg_slice = self.grading_dict.get(deg_tot)
            if deg_slice is None:
                self.grading_dict.update({deg_tot: DegSlice(deg_tot)})
            self.grading_dict[deg_tot].append(vs, idx)
        for (deg, deg_slice) in self.grading_dict.items():
            print(str(deg) + ' '+ str(len(deg_slice.get_vs_list())))
        for (deg, deg_slice) in self.grading_dict.items():
            if not deg_slice.is_complete():
                self.grading_dict.pop(deg)

