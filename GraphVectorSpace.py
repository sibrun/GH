from abc import ABCMeta, abstractmethod
from sage.all import *
import operator
import collections
from tqdm import tqdm
import StoreLoad as SL
import Parallel
import Parameters
import Log
import DisplayInfo

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
    """Vector space interface.

    Abstract class defining the interface for a vector space.

    Attributes:
        properties (VectorSpaceProperties): Vector space properties object, containing information like
            dimension and validity.
    """

    __metaclass__ = ABCMeta

    def __init__(self):
        """Initialize the vector space properties with None."""

        self.properties = VectorSpaceProperties()

    @abstractmethod
    def __eq__(self, other):
        """Comparing two vector spaces.

        Args:
            other (VectorSpace): Vector space to be compared with.

        Returns:
            bool: True if the compared vector spaces are equal, False otherwise.
        """
        pass

    @abstractmethod
    def __str__(self):
        """Unique description of the vector space.

        Returns:
            str: Unique description of the vector space.
        """
        pass

    @abstractmethod
    def get_dimension(self):
        """Dimension of the vector space.

        Returns:
            non-negative int: Dimension of the vector space.
        """
        pass

    @abstractmethod
    def get_ordered_param_dict(self):
        """Get an ordered dictionary of parameters, identifying the vector space.

        Returns:
            Shared.OrderedDict: Ordered dictionary of parameters. Example:
                SH.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops)])
        """
        pass

    @abstractmethod
    def get_work_estimate(self):
        """Estimating the work needed to build the vector space basis.

        Arbitrary unit. Used to schedule the order of building the basis of different vector spaces.

        Returns:
            non-negative int: Estimate the work to build the basis. Arbitrary unit.
        """
        pass

    @abstractmethod
    def build_basis(self):
        """Build the vector space basis."""
        pass

    @abstractmethod
    def update_properties(self):
        """Update the vector space properties."""
        pass

    def get_properties(self):
        """Access the vector space properties.

        Returns:
            VectorSpaceProperties: Vector space properties.
        """
        return self.properties

    def get_sort_dim(self):
        """Dimension for sorting vector spaces.

        Returns:
            non-negative int: Dimension of the vector space if known, constant Parameters.max_sort_value otherwise.
            """
        try:
            sort_dim = self.get_dimension()
        except SL.FileNotFoundError:
            sort_dim = Parameters.max_sort_value
        return sort_dim


class GraphVectorSpace(VectorSpace):
    """Vector space of graphs.

    Abstract class defining the interface for graph vector spaces. It implements the interface vector space and
    provides a method to build the basis.

    Attributes:
        properties (GraphVectorSpaceProperties): Graph vector space properties object, containing information about
            dimension and validity.

    """
    __metaclass__ = ABCMeta

    def __init__(self):
        """Initialize the vector space properties with None."""
        self.properties = GraphVectorSpaceProperties()

    @abstractmethod
    def get_type(self):
        """Returns the type of graphs.

        Returns:
            str: Type of graphs. Example: 'ordinary graphs with evend edges'.
        """
        pass

    @abstractmethod
    def get_basis_file_path(self):
        """Returns the path for the basis file.

        Returns:
            path: Path for the basis file.
        """
        pass

    @abstractmethod
    def get_plot_path(self):
        """Returns the path for the cohomology dimension plot file.

        Returns:
            path: Path for the cohomology dimension plot file.
        """
        pass

    def get_ref_basis_file_path(self):
        """Returns the path for the reference basis file.

        Refers to reference data (if available) for testing.

        Returns:
            path: Path for the reference basis file.
        """
        pass

    @abstractmethod
    def get_partition(self):
        pass

    @abstractmethod
    def is_valid(self):
        """Returns the validity of the parameter combination for the graph vector space.

        Returns:
            bool: True if the graph vector space is valid, False otherwise.
        """
        pass

    @abstractmethod
    def get_generating_graphs(self):
        """Produces a set of graphs whose isomorphism classes span the vector space. (Not necessarily freely!)

        Returns:
            list(sage.Graph): List of sage graphs spanning the vector space.
        """
        pass

    @abstractmethod
    def perm_sign(self, G, p):
        """Returns the sign of the permutation of the edges of graph G, induced by the relabelling by p.

        For G a graph and p a permutation of the edges, returns the sign induced by the relabelling by p.
        Here vertex j becomes vertex p[j] in the new graph.

        Args:
            G (sage.Graph): Sage graph.
            p (list): List of the images of the permutation of the edges of graph G.

        Returns:
            int: Sign of the edge permutation of graph G induced by the relabelling by p.

            """
        pass

    def __str__(self):
        """Unique description of the graph vector space.

        Returns:
            str: Unique description of the grpah vector space.
        """
        return '<%s vector space with parameters: %s>' % (self.get_type(), str(self.get_ordered_param_dict()))

    def graph_to_canon_g6(self, graph):
        """Returns the graph6 string of the canonically labeled graph and the corresponding permutation sign.

        Canonically label the sage Graph graph using the sage method for canonical labelling and respecting the
        partition of the vertices.

        Args:
            graph (sage.Graph): Graph to be canonically labeled.

        Returns:
            tuple(str, int): Tuple containing the graph6 string of the canonically labeled graph and the corresponding
                permutation sign.
        """
        canonG, permDict = graph.canonical_label(partition=self.get_partition(), certificate=True)
        sgn = self.perm_sign(graph, permDict.values())
        return (canonG.graph6_string(), sgn)

    def build_basis(self, progress_bar=False, ignore_existing_files=False, **kwargs):
        """Build the basis of the vector space.

        Create the basis file if the vector space is valid, otherwise skip building a basis. If there exists already
        a basis file rebuild the basis if ignore_existing_file is True, otherwise skip building a basis.
        The basis file contains a list of graph6 strings for canonically labeled graphs building a basis of the vector
        space. The canonical labeling respects the partition of the vertices.

        Args:
            progress_bar (bool, optional): Option to show a progress bar (Default: False).
            ignore_existing_files (bool, optional): Option to ignore existing basis file. Ignore existing file and
                rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
                 (Default: False).
            **kwargs: Accepting further keyword arguments, which have no influence.
        """
        if not self.is_valid():
            # Skip building a basis file if the vector space is not valid.
            return
        if (not ignore_existing_files) and self.exists_basis_file():
            # Skip building a basis file if there exists already one and ignore_existing_file is False.
            return

        generatingList = self.get_generating_graphs()

        if len(generatingList) == 0:
            # Warn if no generating graphs.
            print('Skip building basis, list of generating graphs has zero length for %s' % str(self))
            logger.warn('Skip building basis, list of generating graphs has zero length for %s' % str(self))
            return

        desc = 'Build basis: ' + str(self.get_ordered_param_dict())
        #if not progress_bar:
        print(desc)
        basisSet = set()

        #for G in tqdm(generatingList, desc=desc, disable=(not progress_bar)):
        for G in generatingList:
            # For each graph G in the generating list, add the canonical labeled graph6 representation to the basis set
            # if the graph G doesn't have odd automormphisms.
            if self.get_partition() is None:
                automList = G.automorphism_group().gens()
                canonG = G.canonical_label()
            else:
                # The canonical labelling respects the partition of the vertices.
                automList = G.automorphism_group(partition=self.get_partition()).gens()
                canonG = G.canonical_label(partition=self.get_partition())
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
        """Test whether the graph G has odd automorphisms.

        Args:
            G (sage.Graph): Test whether G has odd automoerphisms.
            automList (list(group generators)): List of generators of the automorphisms group of graph G.

        Returns:
            bool: True if G has odd automorphisms, False otherwise.
            """
        for g in automList:
            if self.perm_sign(G, list(g.tuple())) == -1:
               return True
        return False

    def exists_basis_file(self):
        """Return whether there exists a basis file.

        Returns:
            bool: True if there exists a basis file, False otherwise.
        """
        return os.path.isfile(self.get_basis_file_path())

    def get_dimension(self):
        """Returns the Dimension of the vector space.

        Raises an exception if no basis file found.

        Returns:
            non-negative int: Dimension of the vector space

        Raises:
            StoareLoad.FileNotFoundError: Raised if no basis file found.
        """
        if not self.is_valid():
            return 0
        try:
            header = SL.load_line(self.get_basis_file_path())
            return int(header)
        except SL.FileNotFoundError:
            raise SL.FileNotFoundError("Dimension unknown for %s: No basis file" % str(self))

    def _store_basis_g6(self, basisList):
        """Stores the basis to the basis file.

        The basis file contains a list of graph6 strings for canonically labeled graphs building a basis of the
        vector space.
        The first line of the basis file contains the dimension of the vector space.

        Args:
            basisList (list(str)): List of graph6 strings.
        """
        basisList.insert(0, str(len(basisList)))
        SL.store_string_list(basisList, self.get_basis_file_path())

    def _load_basis_g6(self):
        """Loads the basis from the basis file.

        Raises an exception if no basis file found or if the dimension in the header of the basis file doesn't
        correspond to the dimension of the basis.

        Returns:
            list(graph6 string): List of graph6 strings of canonically labeled graphs building a basis of the vector
                space.

        Raises:
            StoareLoad.FileNotFoundError: Raised if no basis file found.
            ValueError: Raised if dimension in header doesn't correspond to the basis dimension.
        """
        if not self.exists_basis_file():
            raise SL.FileNotFoundError("Cannot load basis, No basis file found for %s: " % str(self))
        basisList = SL.load_string_list(self.get_basis_file_path())
        dim = int(basisList.pop(0))
        if len(basisList) != dim:
            raise ValueError("Basis read from file %s has wrong dimension" % str(self.get_basis_file_path()))
        return basisList

    def get_basis(self, g6=True):
        """Return the basis of the vector space.

        Choose between graph6 and sage graph representation of the basis elements.

        Args:
            g6 (bool): If true a list of graph6 strings is returned. If False a list of sage graphs is returned.

        Returns:
            list(graph6 str / sage.Graph): List of basis elements. As graph6 strings or sage graphs.
        """

        if not self.is_valid():
            logger.warn("Empty basis: %s is not valid" % str(self))
            return []
        basis_g6 = self._load_basis_g6()
        if g6:
            return basis_g6
        else:
            return [Graph(g6) for g6 in basis_g6]

    def get_g6_coordinates_dict(self):
        """Returns a dictionary to translate from the graph6 string of graphs in the basis to their coordinates.

        Returns:
            dict(graph6 str -> int): Dictionary to translate from graph6 string to the coordinate of a basis element.
        """
        return {G6: i for (i, G6) in enumerate(self.get_basis(g6=True))}

    def delete_basis_file(self):
        """Delete the basis file."""
        if os.path.isfile(self.get_basis_file_path()):
            os.remove(self.get_basis_file_path())

    def update_properties(self):
        """Update the graph vector space properties validity and dimension.

        Reading vector space dimension from basis file.

        Raises:
            StoareLoad.FileNotFoundError: Raised if no basis file found.
        """
        self.properties.valid = self.is_valid()
        if not self.properties.valid:
            self.properties.dimension = 0
        else:
            try:
                self.properties.dimension = self.get_dimension()
            except SL.FileNotFoundError:
                pass


class SumVectorSpace(VectorSpace):
    __metaclass__ = ABCMeta

    def __init__(self, vs_list):
        self.vs_list = vs_list
        super(SumVectorSpace, self).__init__()
        self.info_tracker = None

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_ordered_param_range_dict(self):
        pass

    def __eq__(self, other):
        pass

    def get_ordered_param_dict(self):
        pass

    def __str__(self):
        return '<%s vector space with parameters: %s>' % (str(self.get_type()), str(self.get_ordered_param_range_dict()))

    def get_vs_list(self):
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

    '''def __eq__(self, other):
        if len(self.vs_list) != len(other.vs_list):
            return False
        eq_l = 0
        for (vs1, vs2) in itertools.product(self.vs_list, other.vs_list):
            if vs1 == vs2:
                eq_l += 1
        return eq_l == len(self.vs_list)'''

    def sort(self, key='work_estimate'):
        if isinstance(self, DegSlice):
            return
        if key == 'work_estimate':
            self.vs_list.sort(key=operator.methodcaller('get_work_estimate'))
        elif key == 'dim':
            self.vs_list.sort(key=operator.methodcaller('get_sort_dim'))
        else:
            raise ValueError("Invalid sort key. Options: 'work_estimate', 'dim'")

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False, info_tracker=False):
        print(' ')
        print('Build basis of %s' % str(self))
        if n_jobs > 1:
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
        try:
            param_names = self.vs_list[0].get_ordered_param_dict().keys()
            property_names = self.vs_list[0].get_properties().names()
        except IndexError:
            param_names = property_names = []
        parameter_list = param_names + property_names
        self.info_tracker.set_parameter_list(parameter_list)

    def start_tracker(self):
        self.info_tracker = DisplayInfo.InfoTracker(str(self))
        self.set_tracker_parameters()
        vs_info_dict = collections.OrderedDict()
        for vs in self.vs_list:
            vs_info_dict.update({tuple(vs.get_ordered_param_dict().values()): vs.get_properties().list()})
        self.info_tracker.update_data(vs_info_dict)
        self.info_tracker.start()

    def update_tracker(self, vs):
        vs.update_properties()
        message = {tuple(vs.get_ordered_param_dict().values()): vs.get_properties().list()}
        self.info_tracker.get_queue().put(message)

    def stop_tracker(self):
        self.info_tracker.stop()


class DegSlice(SumVectorSpace):
    def __init__(self, vs_list, deg):
        self.deg = deg
        super(DegSlice, self).__init__(vs_list)
        self.start_idx_dict = None

    def __str__(self):
        return '<degree slice with parameters: %s>' % str(self.get_ordered_param_dict())

    def get_type(self):
        pass

    def get_ordered_param_range_dict(self):
        pass

    @abstractmethod
    def __eq__(self, other):
        pass

    @abstractmethod
    def get_ordered_param_dict(self):
        pass

    def build_basis(self, **kwargs):
        super(DegSlice, self).build_basis(**kwargs)
        if not self.is_complete():
            raise ValueError('deg slice %s should be completely built' % str(self))

    def build_start_idx_dict(self):
        self.start_idx_dict = dict()
        start_idx = 0
        for vs in self.vs_list:
            self.start_idx_dict.update({vs: start_idx})
            dim = vs.get_dimension()
            start_idx += dim

    def get_start_idx(self, vector_space):
        if self.start_idx_dict is None:
            self.build_start_idx_dict()
        start_idx = self.start_idx_dict.get(vector_space)
        if start_idx is None:
            raise ValueError('vector_space should refer to a vector space of the degree slice')
        return start_idx

    def is_complete(self):
        if len(self.vs_list) != self.deg + 1:
            return False
        for vs in self.vs_list:
            if vs is None or (vs.is_valid() and not vs.exists_basis_file()):
                return False
        return True
