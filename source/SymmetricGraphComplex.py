""" This file contains auxiliary methods and classes for graph complexes
that are equipped with a symmetric group action.
The typical example is a graph complex with N numbered hairs, with the
symmetric group acting by permutation of the hair labels.
"""

import itertools
from operator import truediv
from sage.all import *
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Log
import Shared
import StoreLoad
import Parameters
import PlotCohomology
from abc import ABCMeta, abstractmethod

logger = Log.logger.getChild('symmetric_graph_complex')


class SymmetricGraphVectorSpace(GraphVectorSpace.GraphVectorSpace):
    """ This abstract class encodes GraphVector spaces with an action of Sn.

    """

    __metaclass__ = ABCMeta

    @abstractmethod
    def get_n(self):
        """ Return n as in S_n, the symmetric group acting on the vector space."""
        pass

    @abstractmethod
    def vertex_permutation_from_permutation(self, p):
        """Returns the permutation on the vertex set of graphs that corresponds to the given permutation p in S_n.
        Hence this method effectively encodes the S_n action.
        : param p: The permutation in S_n, as produced by Permutations(n). In particular, the indices are 1-based.
        : type p: Permutation
        : return: The permutation on vertices as used by GH, i.e., as a list of zero based indices.
        : rtype: list(int)
        """
        pass

    @abstractmethod
    def get_isotypical_projector(self, rep_index):
        """Returns the SymmetricProjectionOperator corresponding to the isotypical component corresponding to
        the rep_index-th irrep (as in Partitions(n))
        """
        pass


class SymmetricProjectionOperator(GraphOperator.GraphOperator):
    """This abstract class encodes the projection operator to an isotypical component of the symmetric group action
        by permuting numbered hairs.
        Warning: The matrix stores not the projector, but projector * c, with c = n!, to have integral matrices.
        The constant c can be obtained by get_normalizing_c().
        The main method to be implemented by the user is

        Deriving
    """

    def __init__(self, domain, rep_index):
        """Initialize the domain and target vector space of the graph operator.

        : param domain: Domain vector space of the operator. This is also the target.
        : type domain: GraphVectorSpace
        : param n: The symmetric group acting should be S_n
        : type n: int
        : param rep_index: The index of the representation in the list produced by Partitions(h).
        : type rep_index: int
        """

        self.rep_index = rep_index
        self.domain = domain
        self.n = domain.get_n()
        n = self.n

        if n <= 0:
            raise ValueError(
                "Error: SymmetricProjectionOperator should not be constructed on vector space with no Sn action.")

        super(SymmetricProjectionOperator, self).__init__(domain, domain)

        # pre-fill in representation and character
        if len(Partitions(n)) <= rep_index:
            raise ValueError(
                "Illegal rep_index: It is larger then the number of irreps.")

        self.rep_partition = Partitions(n)[rep_index]
        self.rep_dim = symmetrica.charvalue(
            self.rep_partition, [1 for j in range(n)])

        self.norm_char_perm = [(symmetrica.charvalue(self.rep_partition, p.cycle_type(
        )), self.domain.vertex_permutation_from_permutation(p)) for p in Permutations(n)]

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding graph operator.
        """
        return domain == target

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        return self.domain.get_dimension()

    def get_type(self):
        return 'projection (irrep=%s)' % str(self.rep_partition)

    def operate_on(self, G):
        # print("operate_on", self)
        # Operates on the graph G with the projection operator
        image = []
        for (c, p) in self.norm_char_perm:
            # c is char value, p is permutation
            # print(p)
            G1 = copy(G)
            sgn = self.domain.perm_sign(G1, p)
            G1.relabel(p, inplace=True)
            image.append((G1, sgn * c))

        return image

    def get_normalizing_c(self):
        """ Returns the normalization constant c, that is, the matrix is c*P, with P the projector.
        """
        return factorial(self.n) / self.rep_dim
        # return factorial(self.n)

    def trace_rank(self):
        """ For projection operators the rank equals the trace. Hence one does not need to compute the rank
        via the default procedure, but can take this faster method as a shortcut. """
        A = self.get_matrix()
        return A.trace() / self.get_normalizing_c()

    def _compute_rank(self, sage=None, linbox=None, rheinfall=None, prime=Parameters.prime):
        # We override _compute_rank so as to avoid expensive rank computation by other means
        if self.is_trivial() or self.get_matrix_entries() == 0:
            return {'exact': 0}
        else:
            return {'exact': self.trace_rank()}


class SymmetricDegSlice(GraphVectorSpace.DegSlice):
    @abstractmethod
    def get_n(self):
        """ Return n as in S_n, the symmetric group acting on the vector space."""
        pass

    @abstractmethod
    def get_isotypical_projector(self, rep_index):
        """Returns the SymmetricProjectionOperator corresponding to the isotypical component corresponding to
        the rep_index-th irrep (as in Partitions(n))
        """
        pass


class SymmetricProjectionOperatorDegSlice(GraphOperator.OperatorMatrix):
    """ Represents the projection operator on a degree slice. 
    """

    def __init__(self, domain, rep_index):
        self.rep_index = rep_index
        self.domain = domain
        self.n = domain.get_n()
        n = self.n

        if n <= 0:
            raise ValueError(
                "Error: SymmetricProjectionOperatorDegSlice should not be constructed on vector space with no Sn action.")

        if len(Partitions(n)) <= rep_index:
            raise ValueError(
                "Illegal rep_index: It is larger then the number of irreps.")

        self.rep_partition = Partitions(n)[rep_index]
        self.rep_dim = symmetrica.charvalue(
            self.rep_partition, [1 for j in range(n)])

        super(SymmetricProjectionOperatorDegSlice,
              self).__init__(domain, domain)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding graph operator.
        """
        return domain == target

    @classmethod
    def generate_op_matrix_list(cls, graded_sum_vs):
        """Return a list of all possible bi operator matrices of this type with domain and target being degree slices
        of the graded sum vector space.

        :param graded_sum_vs: Graded sum vector space composed of degree slices.
        :type graded_sum_vs: GraphVectorSpace.SumVectorSpace
        :return: List of all possible bi operator matrices with domain and target being degree slices of the
            graded sum vector space.
        :rtype: list(BiOperatorMatrix)
        """
        graded_sum_vs_list = graded_sum_vs.get_vs_list()
        bi_op_matrix_list = []
        for (domain, target) in itertools.permutations(graded_sum_vs_list, 2):
            if cls.is_match(domain, target):
                bi_op_matrix_list.append(cls(domain, target))
        return bi_op_matrix_list

    def __str__(self):
        """Return a unique description of the bi operator matrix.

        :return: Unique description of the bi operator matrix.
        :rtype: str
        """
        return '<projection operator (%s): %s>' % (str(self.rep_partition), str(self.domain))

    def is_valid(self):
        return True

    def get_work_estimate(self):
        """Estimate the work needed to build the bi operator matrix by the product of the dimensions of domain and target.

        Used to schedule the order of building the operator matrices.

        :return int: Estimate the work to build the operator matrix.
        :rtype: int
        """

        return self.domain.get_dimension() * self.target.get_dimension()

    def build_matrix(self, ignore_existing_files=False, progress_bar=False, **kwargs):
        if (not ignore_existing_files) and self.exists_matrix_file():
            return
        print(' ')
        print('Build matrix of %s' % str(self))
        shape = (self.domain.get_dimension(), self.target.get_dimension())
        underlying_matrices = self._get_underlying_matrices()
        self._build_underlying_matrices(underlying_matrices, ignore_existing_files=ignore_existing_files,
                                        progress_bar=progress_bar)
        matrix_list = self._get_matrix_list(underlying_matrices)
        self._store_matrix_list(matrix_list, shape)

    def _get_underlying_matrices(self):
        return [vs.get_isotypical_projector(self.rep_index) for vs in self.domain.get_vs_list()]

    def _build_underlying_matrices(self, op_matrix_list, **kwargs):
        for op in op_matrix_list:
            op.build_matrix(**kwargs)

    def _get_matrix_list(self, underlying_matrices):
        matrixList = []
        for op in underlying_matrices:
            if not op.is_valid():
                continue
            domain_start_idx = self.domain.get_start_idx(op.get_domain())
            target_start_idx = self.target.get_start_idx(op.get_target())
            subMatrixList = op.get_matrix_list()
            for (i, j, v) in subMatrixList:
                matrixList.append(
                    (i + domain_start_idx, j + target_start_idx, v))
        matrixList.sort()
        return matrixList


class IsotypicalComponent():
    """Represents an isotypical component in a graph vector space.
    Attributes:
    vs : The underlying vector space
    opP: the projection operator
    rep_index: index of irrep (as in Partitions(n) )
    """

    def __init__(self, vs : SymmetricGraphVectorSpace, rep_index):
        """ vs is the SymmetricGraphVectorSpace, opP the corresponding projecton operator
        rep_index the index of the irrep (as in Partitions(n)) for ehich we take the isotypical component."""
        self.vs : SymmetricGraphVectorSpace = vs
        self.opP = vs.get_isotypical_projector(rep_index)
        self.rep_index = rep_index

    def __eq__(self, other):
        return self.vs == other.vs and self.rep_index == other.rep_index

    def get_dimension(self):
        # for consistency, return dimension of underlying vector space, not of the isotypical component
        return self.vs.get_dimension()

    def get_iso_dimension(self):
        return self.opP.get_matrix_rank()

    def is_valid(self):
        return self.vs.is_valid()

    def __str__(self):
        return "Iso %d of %s" % (self.rep_index, str(self.vs))

    def __hash__(self):
        return hash(str(self))

    def get_ordered_param_dict(self):
        d = self.vs.get_ordered_param_dict().copy()
        d.update({'rep_index': self.rep_index})
        return d
    
    def exists_basis_file(self):
        return self.vs.exists_basis_file()

# class SymmetricSumVectorSpace(GraphVectorSpace.SumVectorSpace):
#     """ Represents the collection of isotypical components."""
#     pass


class SymmetricRestrictedOperatorMatrix(GraphOperator.OperatorMatrix):
    """Represents the restriction of an operator on a symmetric graph complex to an isotypical component.
    For P a projection operator and D a differential the matrix actually stored is
    [D; c(1-P)] to avoid having to compute matrix products.
    """

    def __init__(self, opD, rep_index):
        """ opD represents the differential, opP the projection.
        """
        self.opD = opD
        self.rep_index = rep_index
        self.opP = opD.domain.get_isotypical_projector(rep_index)
        # if opD.domain != opP.domain:
        #     raise ValueError("Domain %s and target %s don't match to build the symmetric composite operator matrix %s"
        #                      % (str(opD.domain), str(self.opP.domain), str(self)))
        super(SymmetricRestrictedOperatorMatrix, self).__init__(
            IsotypicalComponent(opD.domain, rep_index), IsotypicalComponent(opD.target, rep_index))
        # self.is_pseudo_matrix = True

    # @staticmethod
    # def is_match(domain, target):
    #     return True

    def __str__(self):
        return "<Restriction to iso %d of operator %s>" % (self.rep_index, str(self.opD))

    def get_type(self):
        return "Restriction to iso %d of %s" % (self.rep_index, self.opD.get_type())

    def get_work_estimate(self):
        return self.opD.get_work_estimate()

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        """ Builds the matrix by multiplying D and P"""
        # Build matrix of underlying projection operator
        self.opP.build_matrix(ignore_existing_files=ignore_existing_files,
                              skip_if_no_basis=skip_if_no_basis, progress_bar=progress_bar, **kwargs)

        if not self.is_valid():
            return

        if (not ignore_existing_files) and self.exists_matrix_file():
            return

        PM = self.opP.get_matrix_transposed()
        DM = self.opD.get_matrix_transposed()

        myM = PM * DM

        (shapem, shapen) = self.opD.get_matrix_shape()

        myML = [(i, j, myM[i, j])
                for (i, j) in myM.nonzero_positions(copy=False)]

        self._store_matrix_list(myML, (shapen, shapem))

        # if (Dcols != Pcols or Prows != Pcols):
        #     raise ValueError(
        #         "Something wrong: number of columns should match (D is %d x %d, P is %d x %d)!" % (Drows, Dcols, Prows, Pcols))
        # newShape = (Drows+Prows, Dcols)
        # if (Drows != Prows or Prows != Pcols):
        #     raise ValueError(
        #         "Something wrong: number of rows should match (D is %d x %d, P is %d x %d)!" % (Drows, Dcols, Prows, Pcols))
        # newShape = (Drows, Dcols+Pcols)

        # c = self.opP.get_normalizing_c()
        # newList = D_list + [(a, b+Dcols, -v) for (a, b, v)
        #                     in P_list] + [(j, j+Dcols, c) for j in range(Prows)]

        # self._store_matrix_list(newList, newShape)

        # def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        # """ Builds the matrix by stcking [D; 1-P]"""
        # if not self.is_valid():
        #     return

        # # Build matrix of underlying projection operator
        # self.opP.build_matrix(ignore_existing_files=ignore_existing_files,
        #                       skip_if_no_basis=skip_if_no_basis, progress_bar=progress_bar, **kwargs)

        # if (not ignore_existing_files) and self.exists_matrix_file():
        #     return

        # (D_list, Dshape) = self.opD._load_matrix_list()
        # (P_list, Pshape) = self.opP._load_matrix_list()

        # (Drows, Dcols) = Dshape
        # (Prows, Pcols) = Pshape

        # # if (Dcols != Pcols or Prows != Pcols):
        # #     raise ValueError(
        # #         "Something wrong: number of columns should match (D is %d x %d, P is %d x %d)!" % (Drows, Dcols, Prows, Pcols))
        # # newShape = (Drows+Prows, Dcols)
        # if (Drows != Prows or Prows != Pcols):
        #     raise ValueError(
        #         "Something wrong: number of rows should match (D is %d x %d, P is %d x %d)!" % (Drows, Dcols, Prows, Pcols))
        # newShape = (Drows, Dcols+Pcols)

        # c = self.opP.get_normalizing_c()
        # newList = D_list + [(a, b+Dcols, -v) for (a, b, v)
        #                     in P_list] + [(j, j+Dcols, c) for j in range(Prows)]

        # self._store_matrix_list(newList, newShape)

    # def get_matrix_rank(self):
    #     # The rank of the restriction of the operator is not the rank of the matrix actually stored, but the coranks agree.
    #     # Hence the get_matrix_rank method is overridden.
    #     r = super().get_matrix_rank()
    #     fulldim = self.domain.get_dimension()
    #     isodim = self.domain.get_iso_dimension()
    #     print("getrank",r,fulldim,isodim,str(self))
    #     return isodim - (fulldim - r)

    def is_valid(self):
        return self.opD.is_valid()

    def compute_rank(self, sage=None, linbox=None, rheinfall=None, ignore_existing_files=False, skip_if_no_matrix=True):
        print("Compute projector rank "+str(self.opP))
        self.opP.compute_rank(sage, linbox, rheinfall,
                              ignore_existing_files, skip_if_no_matrix)
        print("Done")
        return super().compute_rank(sage, linbox, rheinfall, ignore_existing_files, skip_if_no_matrix)


class SymmetricGraphOperator(GraphOperator.GraphOperator):
    __metaclass__ = ABCMeta

    @abstractmethod
    def restrict_to_isotypical_component(self, rep_index):
        """Returns the SymmetricRestrictedOperatorMatrix that represents the restriction of the operator
        to an isotypical component"""
        pass


class SymmetricBiOperatorMatrix(GraphOperator.BiOperatorMatrix):
    __metaclass__ = ABCMeta

    @abstractmethod
    def restrict_to_isotypical_component(self, rep_index):
        """Returns the SymmetricRestrictedOperatorMatrix that represents the restriction of the operator
        to an isotypical component"""
        pass


class SymmetricDifferential(GraphOperator.Differential):
    """ Represents a differential on a symmetric graph complex, on a per-isotypical-component basis.
    The typical usage is that an ordinary differential holds the SymmetricGraphOperators.
    Then a new SymmetricDifferential is constructed, passing the old differential as constructor paramter.

    Attributes: 
    diff - The (old) differential from which this differential (on isotypical components) is built, 
           by splitting the relevant operators into their restrictions on istypical components.
    """

    def __init__(self, diff):
        """ Initializes the RestrictedContractEdgesD-differential from a ContractEdgesD object.
        Before construction, cohomology for ContractEdgesD should be available, since we will add only those
        operators that are necessary for computing nonzero cohomology."""
        self.diff = diff
        (vsList, opList) = SymmetricDifferential.split_isotypical_components(diff)
        super(SymmetricDifferential, self).__init__(
            GraphVectorSpace.SumVectorSpace(vsList), opList)

    def refine_cohom_dim_dict(self, dict):
        """Refines the given dictionary of cohomology dimensions (vectorspace->int) by providing info on the splitting into
        Sn irreducible components. 
        Returns a new dictionary, the old is unchanged. """
        newdict = dict.copy()
        mydict = self._get_cohomology_dim_dict()
        # print(mydict)
        for vs, dim in dict.items():
            # look for available refinements
            refines = []
            for iso, val in mydict.items():
                if iso.vs == vs:
                    # found match
                    refines.append("%d$s_{%s}$" %
                                   (val, str(iso.opP.rep_partition)))
            if len(refines) > 0:
                newdim = str(dim) + " (" + "+".join(refines) + ")"
                newdict[vs] = newdim
        return newdict

    def _get_cohomology_dim_dict(self):
        # need tooverride this to correct for dimensions
        d = super()._get_cohomology_dim_dict()
        for vs, dim in d.items():
            print(".... "+str(vs))
            print(vs.get_dimension(), vs.get_iso_dimension())
            d[vs] = d[vs] - vs.get_dimension() + vs.get_iso_dimension()
        return d

    def plot_cohomology_dim(self, to_html=False, to_csv=False, x_plots=2):
        self.plot_refined_cohomology_dim(to_html, to_csv, x_plots)

    def plot_refined_cohomology_dim(self, to_html=False, to_csv=False, x_plots=2):
        """Plot the cohomology dimensions, including data on isotypical decompositions.

        Plot the cohomology dimensions as plot and/or table associated with the differential.

        :param to_html: Option to generate a html file with a table of the cohomology dimensions (Dafault: False).
        :type to_html: bool
        :param to_csv: Option to generate a csv file with a table of the cohomology dimensions (default: False).
        :type to_csv: bool
        :param x_plots: Number of plots on the x-axis (Default: 2).
        :type x_plots: int
        """
        print(' ')
        print('Plot cohomology dimensions of the associated graph complex of ' + str(self.diff))
        logger.warning(
            'Plot cohomology dimensions of the associated graph complex of ' + str(self.diff))
        dim_dict = self.diff._get_cohomology_dim_dict()
        dim_dict_refined = self.refine_cohom_dim_dict(dim_dict)

        # look up parameter values
        dim_dict_refined2 = dict()
        for vs in self.diff.sum_vector_space.get_vs_list():
            dim_dict_refined2.update(
                {vs.get_ordered_param_dict().get_value_tuple(): dim_dict_refined.get(vs)})

        plot_path = self.get_cohomology_plot_path()
        # print(plot_path)
        parameter_order = self.diff.get_cohomology_plot_parameter_order()
        ordered_param_range_dict = self.diff.get_ordered_cohomology_param_range_dict()
        PlotCohomology.plot_array(dim_dict_refined2, ordered_param_range_dict, plot_path, to_html=to_html, to_csv=to_csv,
                                  x_plots=x_plots, parameter_order=parameter_order)

    @classmethod
    def split_isotypical_components(cls, diff):
        """Creates an operator list from differential diff by restricting to the isotypical components.
        Only those operators are added that are involved in nonzero cohomology entries, with action of Sn with n>=2, to avoid duplicate computations.
        Also, if the irrep has too high dimension, the isotypical component is not considered either.

        It is assumed that matrix ranks of the original operators have been computed before calling this method.

        From the lists returned an SymmetricDifferential will usually be constructed, with the operator collection build from SymmetricRestrictedOperatorMatrix
        objects.

        Returns:
        A pair (vsList, opList) with vsList the list of isotypical components, opList the list of (restricted) operators
        """
        opD_list = []
        for (opD, opDD) in itertools.permutations(diff.op_matrix_list, 2):
            if opD.get_domain() == opDD.get_target():
                n = opD.get_domain().get_n()
                if n >= 2:
                    dim = GraphOperator.Differential.cohomology_dim(opD, opDD)
                    #dim = 1
                    if type(dim) == int and dim > 0:
                        # We have found non-zero cohomology.

                        # gather restricted differentials
                        for rep_ind in range(len(Partitions(opD.domain.get_n()))):
                            opA = opD.restrict_to_isotypical_component(rep_ind)
                            if opA.opP.rep_dim > dim:
                                continue
                            if not (opA in opD_list):
                                opD_list.append(opA)

                            opA = opDD.restrict_to_isotypical_component(
                                rep_ind)
                            if not (opA in opD_list):
                                opD_list.append(opA)
        vsList = [op.domain for op in opD_list]
        return (vsList, opD_list)


# def getCohomDimP(op1, op2, opP, rep_ind):
#     """
#     Computes the cohomology of the isotypical component corresponding to the projector P of the complex.
#     :param op1: The operator corresponding to the differential D.
#     : type op1: GraphOperator
#     :param op2: The operator corresponding to the differential DD, the cohomology is kerD/imDD.
#     : type op2: GraphOperator
#     :param opP: The projection operator corresponding to the isotypical component (i.e. im P).
#     : type op1: GraphOperator
#     : return : The
#     """
#     # tt = CHairyGraphComplex.ContractEdgesGO.generate_operator(
#     #     n_vertices, n_loops, n_hairs, even_edges)
#     # tu = CHairyGraphComplex.ContractEdgesGO.generate_operator(
#     #     n_vertices+1, n_loops, n_hairs, even_edges)
#     # symmp1 = CHairyGraphComplex.SymmProjector.generate_operator(
#     #     n_vertices, n_loops, n_hairs, even_edges, rep_ind)

#     D1 = op1.get_matrix()
#     D2 = op2.get_matrix()
#     # C = D2*D1
#     opP.build_matrix(ignore_existing_files=True)
#     P1 = opP.get_matrix()
#     print("matrices loaded")

#     D1P = D1*P1
#     D2P = P1*D2
#     print("computing ranks....")

#     # diff = D1*D2
#     # diffs = sum(abs(c) for cc in diff.columns() for c in cc)
#     # print(diffs)

#     isocomp_dim = P1.rank()  # todo: take trace here
#     r1 = D1P.rank()
#     r2 = D2P.rank()
#     # print(isocomp_dim, r1, r2)
#     cohomdim = isocomp_dim - r1-r2
#     n = opP.domain.get_n()
#     part = Partitions(n)[rep_ind]
#     rep_dim = symmetrica.charvalue(part, [1 for j in range(n)])
#     # if cohomdim > 0:
#     #     print("Cohomology found:  ", "even_edges" if even_edges else "odd_edges", ", h=", n_hairs, ", l=", n_loops, ", vertices=", n_vertices,
#     #           " (degree ", n_vertices+n_loops -
#     #           1, "), partition=", part,  ", invpartition=", part.conjugate(),
#     #           ", multiplicity=", cohomdim/rep_dim, ", cohomdim=", cohomdim)
#     # return isocomp_dim - r1-r2
#     return cohomdim / rep_dim


# def getCohomDimPAll(gvs):
#     for h in gvs.h_range:
#         for l in gvs.l_range:
#             for v in gvs.v_range:
#                 D1 = CHairyGraphComplex.ContractEdgesGO.generate_operator(
#                     v, l, h, gvs.even_edges)
#                 D2 = CHairyGraphComplex.ContractEdgesGO.generate_operator(
#                     v+1, l, h, gvs.even_edges)
#                 try:
#                     d = CHairyGraphComplex.CHairyGraphVS(
#                         v, l, h, gvs.even_edges).get_dimension()
#                     r1 = D1.get_matrix_rank()
#                     r2 = D2.get_matrix_rank()
#                     if d-r1-r2 > 0:
#                         for rep_ind in range(len(Partitions(h))):
#                             getCohomDimP(v, l, h, gvs.even_edges, rep_ind)
#                 except:
#                     pass
