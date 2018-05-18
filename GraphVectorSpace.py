"""Module providing abstract classes for vector spaces, graph vector spaces, direct sum of vector spaces and
degree slices for bigraded vector spaces."""


__all__ = ['VectorSpaceProperties', 'GraphVectorSpaceProperties', 'VectorSpace', 'GraphVectorSpace', 'SumVectorSpace',
           'DegSlice']

from abc import ABCMeta, abstractmethod
from sage.all import *
import operator
import collections
from tqdm import tqdm
import StoreLoad
import Parallel
import Parameters
import Log
import DisplayInfo


logger = Log.logger.getChild('graph_vector_space')


class VectorSpaceProperties(object):
    """Properties of a vector space.

    Attributes:
        dimension (non-negative int): Dimension of the vector space.
    """
    def __init__(self):
        """Initialize the vector space properties with None."""
        self.dimension = None

    @classmethod
    def names(cls):
        """Returns a list of the names of the vector space properties.

        :return: list(str): Names of the vector space properties.
        """
        return ['dimension']

    @staticmethod
    def sort_variables():
        """Returns a list of the vector space properties used as sort keys for vector spaces.

        :return: list(str): Names of the vector space properties, used as sort keys for vector spaces matrices.
        """
        return VectorSpaceProperties.names()

    def list(self):
        """Returns a list of the vector space properties.

        :return: list: Vector space properties.
        """
        return [self.dimension]


class GraphVectorSpaceProperties(VectorSpaceProperties):
    """Properties of a graph vector space.

    Attributes:
        valid (bool): Validity of the parameter combination of the graph vector space.

        dimension (non-negative int): Dimension of the graph vector space.
    """
    def __init__(self):
        """Initialize the graph vector space properties with None."""
        self.valid = None
        self.dimension = None

    @classmethod
    def names(cls):
        """Returns a list of names of the graph vector space properties.

        :return: list(str): Names of the graph vector space properties.
        """
        return ['valid', 'dimension']

    def list(self):
        """Returns a list of the graph vector space properties.

        :return: list: Graph vector space properties.
        """
        return [self.valid, self.dimension]


class VectorSpace(object):
    """Vector space interface.

    Abstract class defining the interface for a vector space.

    Attributes:
        properties (VectorSpaceProperties): Vector space properties, containing the dimension.
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        """Initialize the vector space properties with None."""
        self.properties = VectorSpaceProperties()

    @abstractmethod
    def __eq__(self, other):
        """Comparing two vector spaces.

        :param other: VectorSpace: Vector space to be compared with.
        :return: bool: True if the compared vector spaces are equal.
        """
        pass

    @abstractmethod
    def __str__(self):
        """Unique description of the vector space.

        :return: str: Unique description of the vector space.
        """
        pass

    @abstractmethod
    def get_dimension(self):
        """Returns the dimension of the vector space.

        :return non-negative int: Dimension of the vector space.
        """
        pass

    @abstractmethod
    def get_ordered_param_dict(self):
        """Returns an ordered dictionary of parameters, identifying the vector space.

        :return: Shared.OrderedDict: Ordered dictionary of parameters. Example:
                SH.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])
        """
        pass

    @abstractmethod
    def get_work_estimate(self):
        """Estimates the work needed to build the vector space basis.

        Arbitrary units. Used to schedule the order of building the basis of different vector spaces.

        :return: non-negative int: Estimate the work to build the basis. Arbitrary units.
        """
        pass

    @abstractmethod
    def build_basis(self, progress_bar=False, ignore_existing_files=False, **kwargs):
        """Build the vector space basis.

        :param progress_bar: bool, optional: Option to show a progress bar (Default: False).
        :param ignore_existing_files: bool, optional: Option to ignore existing basis file. Ignore existing file and
                rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
                 (Default: False).
        :param kwargs: Accepting further keyword arguments.
        """
        pass

    @abstractmethod
    def update_properties(self):
        """Update the vector space properties."""
        pass

    def get_properties(self):
        """Returns the vector space properties.

        :return: VectorSpaceProperties: Vector space properties.
        """
        return self.properties

    def get_sort_dim(self):
        """Dimension for sorting vector spaces.

        :return:non-negative int: Dimension of the vector space if known, constant Parameters.max_sort_value otherwise.
        """
        try:
            sort_dim = self.get_dimension()
        except StoreLoad.FileNotFoundError:
            sort_dim = Parameters.max_sort_value
        return sort_dim


class GraphVectorSpace(VectorSpace):
    """Vector space of graphs.

    Abstract class defining the interface for graph vector spaces. It implements the interface vector space and
    provides a method to build the basis.

    Attributes:
        properties (GraphVectorSpaceProperties): Graph vector space properties, containing information about
            dimension and validity.

    """
    __metaclass__ = ABCMeta

    def __init__(self):
        """Initialize the vector space properties with None."""
        self.properties = GraphVectorSpaceProperties()

    @abstractmethod
    def get_type(self):
        """Returns a unique description of the graph type.

        :return: str: Type of graphs. Example: 'ordinary graphs with even edges'.
        """
        pass

    @abstractmethod
    def get_basis_file_path(self):
        """Returns the path to the basis file.

        :return: path: Path to the basis file.
        """
        pass

    def get_ref_basis_file_path(self):
        """Returns the path to the reference basis file.

        Refers to reference data (if available) for testing.

        :return: path: Path to the reference basis file.
        """
        pass

    @abstractmethod
    def get_partition(self):
        """Returns the partition of the vertices in different colours.

        The partition of vertices is respected for the canonical labelling as well as for the automorphism group of
        graphs.

        :return: list(list(non-negative int): List of lists. Partition of the vertices in different colours. Example:
                [list(range(0, self.n_vertices)), list(range(self.n_vertices, self.n_vertices + self.n_hairs))]
        """
        pass

    @abstractmethod
    def is_valid(self):
        """Returns the validity of the parameter combination for the graph vector space.

        :return: bool: True if the graph vector space is valid.
        """
        pass

    @abstractmethod
    def get_generating_graphs(self):
        """Returns a set of graphs whose isomorphism classes span the graph vector space. (Not necessarily freely!)

        :return: list(sage.Graph): List of sage graphs spanning the graph vector space.
        """
        pass

    @abstractmethod
    def perm_sign(self, G, p):
        """Returns the sign of the permutation of the edges of the graph G, induced by the relabelling by
        the permutation p.

        For G a graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
        Here vertex j becomes vertex p[j] in the new graph.

        :param G: sage.Graph: Sage graph.
        :param p: list(non-negative int): List of the images of the permutation of the edges of graph G.
        :return: int: Sign of the edge permutation of graph G induced by the relabelling by p.
        """
        pass

    def __str__(self):
        """Unique description of the graph vector space.

        :return: str: Unique description of the graph vector space.
        """
        return '<%s vector space with parameters: %s>' % (self.get_type(), str(self.get_ordered_param_dict()))

    def graph_to_canon_g6(self, graph):
        """Returns the graph6 string of the canonically labeled graph and the corresponding permutation sign.

        Labels the sage Graph graph canonically using the sage method for canonical labelling and respecting the
        partition of the vertices.

        :param graph: sage.Graph: Graph to be canonically labeled.
        :return: tuple(str, int): Tuple containing the graph6 string of the canonically labeled graph and the
            corresponding permutation sign.
        """
        canonG, perm_dict = graph.canonical_label(partition=self.get_partition(), certificate=True)
        sgn = self.perm_sign(graph, perm_dict.values())
        return (canonG.graph6_string(), sgn)

    def build_basis(self, progress_bar=False, ignore_existing_files=False, **kwargs):
        """Build the basis of the vector space.

        Creates the basis file if the vector space is valid, otherwise skips building a basis. If there exists already
        a basis file rebuild the basis if ignore_existing_file is True, otherwise skip building a basis.
        The basis file contains a list of graph6 strings for canonically labeled graphs building a basis of the vector
        space. The canonical labeling respects the partition of the vertices.

        :param progress_bar: bool, optional: Option to show a progress bar (Default: False).
        :param ignore_existing_files: bool, optional: Option to ignore existing basis file. Ignore existing file and
                rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
                (Default: False).
        :param kwargs: Accepting further keyword arguments, which have no influence.
        """
        if not self.is_valid():
            # Skip building a basis file if the vector space is not valid.
            return
        if (not ignore_existing_files) and self.exists_basis_file():
            # Skip building a basis file if there exists already one and ignore_existing_file is False.
            return

        generatingList = self.get_generating_graphs()

        desc = 'Build basis: ' + str(self.get_ordered_param_dict())
        #if not progress_bar:
        print(desc)
        basis_set = set()
        #for G in tqdm(generatingList, desc=desc, disable=(not progress_bar)):
        for G in generatingList:
            # For each graph G in the generating list, add the canonical labeled graph6 representation to the basis set
            # if the graph G doesn't have odd automormphisms.
            if self.get_partition() is None:
                autom_list = G.automorphism_group().gens()
                canonG = G.canonical_label()
            else:
                # The canonical labelling respects the partition of the vertices.
                autom_list = G.automorphism_group(partition=self.get_partition()).gens()
                canonG = G.canonical_label(partition=self.get_partition())
            canon6 = canonG.graph6_string()
            if not canon6 in basis_set:
                if not self._has_odd_automorphisms(G, autom_list):
                    basis_set.add(canon6)

        self._store_basis_g6(list(basis_set))

    def _has_odd_automorphisms(self, G, automList):
        """Returns whether the graph G has odd automorphisms.

        :param G: sage.Graph: Test whether G has odd automoerphisms.
        :param automList: list(group generators): List of generators of the automorphisms group of graph G.
        :return: bool: True if G has odd automorphisms.
        """
        for g in automList:
            if self.perm_sign(G, list(g.tuple())) == -1:
               return True
        return False

    def exists_basis_file(self):
        """Return whether there exists a basis file.

        :return: bool: True if there exists a basis file, False otherwise.
        """
        return os.path.isfile(self.get_basis_file_path())

    def get_dimension(self):
        """Returns the Dimension of the vector space.

        :return: non-negative int: Dimension of the vector space

        :raise StoreLoad.FileNotFoundError: Raised if no basis file found.
        """
        if not self.is_valid():
            return 0
        try:
            header = StoreLoad.load_line(self.get_basis_file_path())
            return int(header)
        except StoreLoad.FileNotFoundError:
            raise StoreLoad.FileNotFoundError("Dimension unknown for %s: No basis file" % str(self))

    def _store_basis_g6(self, basis_list):
        """Stores the basis to the basis file.

        The basis file contains a list of graph6 strings for canonically labeled graphs building a basis of the
        vector space.
        The first line of the basis file contains the dimension of the vector space.

        :param basis_list: list(str): List of graph6 strings representing the vector space basis.
        """
        basis_list.insert(0, str(len(basis_list)))
        StoreLoad.store_string_list(basis_list, self.get_basis_file_path())

    def _load_basis_g6(self):
        """Loads the basis from the basis file.

        Raises an exception if no basis file found or if the dimension in the header of the basis file doesn't
        correspond to the dimension of the basis.

        :return: list(graph6 string): List of graph6 strings of canonically labeled graphs building a basis of the
            vector space.
        :raise StoreLoad.FileNotFoundError: Raised if no basis file found.
        :raise ValueError: Raised if dimension in header doesn't correspond to the basis dimension.
        """
        if not self.exists_basis_file():
            raise StoreLoad.FileNotFoundError("Cannot load basis, No basis file found for %s: " % str(self))
        basis_list = StoreLoad.load_string_list(self.get_basis_file_path())
        dim = int(basis_list.pop(0))
        if len(basis_list) != dim:
            raise ValueError("Basis read from file %s has wrong dimension" % str(self.get_basis_file_path()))
        return basis_list

    def get_basis_g6(self):
        """Returns the basis of the vector space as list of graph6 strings.

        :return: list(graph6 str): List of graph6 strings representing the basis elements.
        """
        if not self.is_valid():
            # Return empty list if graph vector space is not valid.
            logger.warn("Empty basis: %s is not valid" % str(self))
            return []
        return self._load_basis_g6()

    def get_basis(self):
        """Returns the basis of the vector space as list of sage graphs.

        :return: list(sage.Graph): List of sage graphs representing the basis elements.
        """
        return map(Graph, self.get_basis_g6())

    def get_g6_coordinates_dict(self):
        """Returns a dictionary to translate from the graph6 string of graphs in the basis to their coordinates.

        :return: dict(graph6 str -> int): Dictionary to translate from graph6 string to the coordinate of a basis element.
        """
        return {G6: i for (i, G6) in enumerate(self.get_basis_g6())}

    def delete_basis_file(self):
        """Delete the basis file."""
        if os.path.isfile(self.get_basis_file_path()):
            os.remove(self.get_basis_file_path())

    def update_properties(self):
        """Update the graph vector space properties validity and dimension.

        Reading vector the space dimension from basis file.

        :raise StoreLoad.FileNotFoundError: Raised if no basis file found.
        """
        self.properties.valid = self.is_valid()
        if not self.properties.valid:
            self.properties.dimension = 0
        else:
            try:
                self.properties.dimension = self.get_dimension()
            except StoreLoad.FileNotFoundError:
                pass


class SumVectorSpace(VectorSpace):
    """Direct sum of vector spaces.

    Abstract class. Implements the interface vector space.

    Attributes:
        vs_list (list(VectorSpace)): List of sub vector spaces.

        info_tracker (DisplayInfo.InfoTracker): Tracker for information about the vector spaces in vs_list.
            Tracker is only active if the basis of different vector spaces are not built in parallel.

        properties (VectorSpaceProperties): Vector space properties object, containing information about the
            dimension.
    """
    __metaclass__ = ABCMeta

    def __init__(self, vs_list):
        """Initialize with a list of vector spaces.

        :param vs_list: list(VectorSpace): List of vector spaces to initialize the sum vector space.
        """
        self.vs_list = vs_list
        self.info_tracker = None
        super(SumVectorSpace, self).__init__()

    def __eq__(self, other):
        pass

    @abstractmethod
    def get_type(self):
        """Returns a unique description of the type of graphs.

        :return str: Type of graphs. Example: 'ordinary graphs with even edges'.
        """
        pass

    @abstractmethod
    def get_ordered_param_range_dict(self):
        """Returns an ordered dictionary of parameter ranges, for the sub vector spaces of the sum vector space.

         :return Shared.OrderedDict: Ordered dictionary of parameter ranges. Example:
                 Shared.OrderedDict([('vertices', self.v_range), ('loops', self.l_range)])
         """
        pass

    def get_ordered_param_dict(self):
        """Returns an ordered dictionary of parameters, identifying the vector space.

         :return Shared.OrderedDict: Ordered dictionary of parameters. Example:
                 Shared.OrderedDict([('deg', self.deg)])
         """
        pass

    def __str__(self):
        """Unique description of the vector space.

        :return str: Unique description of the vector space.
        """
        return '<%s vector space with parameters: %s>' % (str(self.get_type()), str(self.get_ordered_param_range_dict()))

    def get_vs_list(self):
        """Returns the list of sub vector spaces.

        :return list(VectorSpace): List of sub vector spaces.
        """
        return self.vs_list

    def get_work_estimate(self):
        """Estimates the work needed to build the vector space basis.

        Arbitrary units. Used to schedule the order of building the basis of different vector spaces.

        :return non-negative int: Estimate the work to build the basis. Arbitrary units.
        """
        work_estimate = 0
        for vs in self.vs_list:
            work_estimate += vs.get_work_estimate()
        return work_estimate

    def get_dimension(self):
        """Returns the dimension of the sum vector space.

        :return non-negative int: Dimension of the sum vector space.
        """
        dim = 0
        for vs in self.vs_list:
            dim += vs.get_dimension()
        return dim

    def update_properties(self):
        """Update the graph vector space property dimension.

        Reading vector space dimension from basis files.

        :raise StoreLoad.FileNotFoundError: Raised if a basis file is not found.
        """
        try:
            self.properties.dimension = self.get_dimension()
        except StoreLoad.FileNotFoundError:
            pass

    def contains(self, vector_space):
        """Determines whether the vector space vector_space is a sub vector space of the sum vector space.

        :param vector_space: VectorSpace: Test whether this vector space is contained in the sum vector space.
        :return: bool: True if the vector space vector_space is a sub vector space of the sum vector space.
        """

        for vs in self.vs_list:
            if vs == vector_space:
                return True
        return False

    def sort(self, key='work_estimate'):
        """Sort the sub vector spaces to schedule building the basis.

        :param key: str, optional: Choose between sort key 'work_estimate' and 'dim' for dimension (Default: 'work_estimate')
        :raise ValueError: Raises exception if sort key is neither 'work_estimate' nor 'dim'
        """
        if isinstance(self, DegSlice):
            return
        if key == 'work_estimate':
            self.vs_list.sort(key=operator.methodcaller('get_work_estimate'))
        elif key == 'dim':
            self.vs_list.sort(key=operator.methodcaller('get_sort_dim'))
        else:
            raise ValueError("Invalid sort key. Options: 'work_estimate', 'dim'")

    def build_basis(self, ignore_existing_files=False, n_jobs=1, progress_bar=False, info_tracker=False):
        """Build the basis of the sub vector spaces.

        Call build_basis for each sub vector space of the sum vector space.

        :param ignore_existing_files: bool, optional: Option to ignore existing basis files. Ignore existing files and
                rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
                 (Default: False).
        :param n_jobs: positive int, optional: Option to compute the basis of the different sub vector spaces in parallel using
                n_jobs parallel processes (Default: 1).
        :param progress_bar: bool, optional: Option to show a progress bar (Default: False). Only active if the basis of
                different sub vector spaces ar not built in parallel.
        :param info_tracker: bool, optional: Option to plot information about the sub vector spaces in a web page.
                Only active if basis not built in parallel processes (Default: False).
        """
        print(' ')
        print('Build basis of %s' % str(self))
        if n_jobs > 1:
            # If mor than 1 process progress bar and info tracker are not activated.
            progress_bar = False
            info_tracker = False
        if info_tracker:
            self.start_tracker()
        self.sort()
        Parallel.parallel(self._build_single_basis, self.vs_list, n_jobs=n_jobs, progress_bar=progress_bar,
                        ignore_existing_files=ignore_existing_files, info_tracker=info_tracker)
        if info_tracker:
            self.stop_tracker()

    def _build_single_basis(self, vs, progress_bar=False, ignore_existing_files=True, info_tracker=False):
        info = info_tracker if Parameters.second_info else False
        vs.build_basis(progress_bar=progress_bar, ignore_existing_files=ignore_existing_files, info_tracker=info)
        if info_tracker:
            self.update_tracker(vs)

    def set_tracker_parameters(self):
        """Initialize the info tracker by setting its parameters.

        Set the names of the variables to be displayed.
        """
        try:
            param_names = self.vs_list[0].get_ordered_param_dict().keys()
            property_names = self.vs_list[0].get_properties().names()
        except IndexError:
            param_names = property_names = []
        parameter_list = param_names + property_names
        self.info_tracker.set_parameter_list(parameter_list)

    def start_tracker(self):
        """Start the info tracker.

        Track information about the properties of the sub vector spaces and display it in a web page.
        """
        self.info_tracker = DisplayInfo.InfoTracker(str(self))
        self.set_tracker_parameters()
        vs_info_dict = collections.OrderedDict()
        for vs in self.vs_list:
            vs_info_dict.update({tuple(vs.get_ordered_param_dict().values()): vs.get_properties().list()})
        self.info_tracker.update_data(vs_info_dict)
        self.info_tracker.start()

    def update_tracker(self, vector_space):
        """Update info tracker for the vector space vector_space.

        :param vector_space: VectorSpace: Vector space for which to update the properties and message it to the info tracker.
        """
        vector_space.update_properties()
        message = {tuple(vector_space.get_ordered_param_dict().values()): vector_space.get_properties().list()}
        self.info_tracker.get_queue().put(message)

    def stop_tracker(self):
        """Stop tracking informations about the sub vector spaces."""
        self.info_tracker.stop()


class DegSlice(SumVectorSpace):
    """Special type of a direct sum of vector spaces, to be used as degree slices in bicomplexes, i.e. a bigraded
    vector space is built as a direct sum of degree slices.

    Abstract class. Implementing SumVectorSpace.

    Attributes:
        deg (int): Degree of the degree slice.

        vs_list (list(VectorSpace)): List of sub vector spaces.

        info_tracker (DisplayInfo.InfoTracker): Tracker for information about the vector spaces in vs_list.
            Tracker is only active if the basis of different vector spaces are not built in parallel.

        properties (VectorSpaceProperties): Vector space properties object, containing information about the
            dimension.
    """

    def __init__(self, vs_list, deg):
        """Initialize with a list of vector spaces and a degree.

        :param vs_list: list(VectorSpace): List of vector spaces to initialize the sum vector space.
        :param deg: int: Degree of the degree slice.
        """
        self.deg = deg
        super(DegSlice, self).__init__(vs_list)
        self.start_idx_dict = None

    def __str__(self):
        """Unique description of the vector space.

        :return str: Unique description of the vector space.
        """
        return '<degree slice with parameters: %s>' % str(self.get_ordered_param_dict())

    def get_type(self):
        pass

    def get_ordered_param_range_dict(self):
        pass

    @abstractmethod
    def __eq__(self, other):
        """Comparing two degree slices.

        :param other: DegSlice: Degree slice to be compared with.
        :return: bool: True if the compared degree slices are equal.
        """
        pass

    @abstractmethod
    def get_ordered_param_dict(self):
        """Returns an ordered dictionary of parameters, identifying the degree slice.

         :return Shared.OrderedDict: Ordered dictionary of parameters. Example:
                 Shared.OrderedDict([('deg', self.deg)])
         """
        pass

    def build_basis(self, **kwargs):
        """Build the basis of the sub vector spaces of the degree slice.

        :param kwargs: Forward keword arguments to the build basis method of the SumVectorSpace.

        :raise ValueError: If the basis of the degree slice is not completely built, i.e. not for all valid sub vector
                spaces there exists a basis file.
        """
        super(DegSlice, self).build_basis(**kwargs)
        if not self.is_complete():
            raise ValueError('Degree slice %s should be completely built' % str(self))

    def build_start_idx_dict(self):
        """Builds a dictionary of start indices of the sub vector spaces.

        The dictionary contains the coordinate of the first basis vector of each sub vector space as basis vector of the
        degree slice.
        """
        self.start_idx_dict = dict()
        start_idx = 0
        for vs in self.vs_list:
            self.start_idx_dict.update({vs: start_idx})
            dim = vs.get_dimension()
            start_idx += dim

    def get_start_idx(self, vector_space):
        """Returns the start index of the sub vector space vector_space.

        Returns the coordinate of the first basis vector of the sub vector space as basis vector of the degree slice.

        :param vector_space: VectorSpace: Sub vector space of which to get the start index.
        :return: non-negative int: Start index of the sub vector space vector_space.
        :raise ValueError: If vector_space is not a sub vector space of the degree slice.
        """
        if self.start_idx_dict is None:
            self.build_start_idx_dict()
        start_idx = self.start_idx_dict.get(vector_space)
        if start_idx is None:
            raise ValueError('vector_space should refer to a sub vector space of the degree slice')
        return start_idx

    def is_complete(self):
        """Test whether the degree slice is complete.

        For all valid sub vector spaces of the degree slice a basis file should exist.

        :return bool: True if the degree slice is complete, False otherwise.
        """
        if len(self.vs_list) != self.deg + 1:
            return False
        for vs in self.vs_list:
            if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
                return False
        return True
