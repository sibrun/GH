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
import Shared
import StoreLoad
import Parameters


class SymmetricGraphVectorSpace(GraphVectorSpace.GraphVectorSpace):
    """ This abstract class encodes GraphVector spaces with an action of Sn.
    
    """

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





class SymmetricProjectionOperator(GraphOperator.GraphOperator):
    """This abstract class encodes the projection operator to an isotypical component of the symmetric group action
        by permuting numbered hairs.
        Warning: The matrix stores not the projector, but projector * n! / rep_dimension, to have integral matrices.

        The main method to be implemented by the user is 

        Deriving 
    """



    def __init__(self, domain, rep_index):
        """Initialize the domain and target vector space of the contract edges graph operator.

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

        super(SymmetricProjectionOperator, self).__init__(domain, domain)

        # pre-fill in representation and character
        self.rep_partition = Partitions(n)[rep_index]
        self.norm_char_perm = [(symmetrica.charvalue(self.rep_partition, p.cycle_type(
            )), self.domain.vertex_permutation_from_permutation(p)) for p in Permutations(n)]

        # print(self.norm_char_perm)

    @staticmethod
    def is_match(domain, target):
        """Check whether domain and target match to generate a corresponding graph operator.
        """
        return domain == target

    def get_work_estimate(self):
        # Returns as work estimate: domain.n_edges * log(target dimension, 2)
        return self.domain.get_dimension()

    def get_type(self):
        return 'projection operator'

    def operate_on(self, G):
        # Operates on the graph G with the projection operator
        image = []
        for (c, p) in self.norm_char_perm:
            # c is char value, p is permutation
            G1 = copy(G)
            sgn = self.domain.perm_sign(G1, p)
            G1.relabel(p, inplace=True)
            image.append((G1, sgn * c))

        return image



class SymmetricCompositeOperatorMatrix(GraphOperator.OperatorMatrix):
    """Represents the restriction of an operator on a symmetric graph complex to an isotypical component. 
    More precisely, for P a projection operator and D a differential the matrix stored is 
    [D; 1-P].
    """

    def __init__(self, opD, opP):
        """ opD represents the differential, opP the projection. 
        """
        self.opD = opD
        self.opP = opP 
        if opD.domain != opP.domain:
            raise ValueError("Domain %s and target %s don't match to build the symmetric composite operator matrix %s"
                             % (str(opD.domain), str(opP.domain), str(self)))

        super(SymmetricCompositeOperatorMatrix, self).__init__(opD.domain, opD.target)

    # @staticmethod
    # def is_match(domain, target):
    #     return True

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        """ Builds the matrix by stcking [D; 1-P]"""
        if not self.is_valid():
            return
        if (not ignore_existing_files) and self.exists_matrix_file():
            return
        

        (D_list, Dshape) = self.opD._load_matrix_list()
        (P_list, Pshape) = self.opP._load_matrix_list()

        (Drows, Dcols) = Dshape
        (Prows, Pcols) = Pshape

        if (Dcols != Pcols or Prows != Pcols):
            raise ValueError("Something wrong: number of columns should match!")
        newShape = (Drows+Prows, Dcols)

        newList = D_list + [(a+Drows,b,-c) for (a,b,c) in P_list] + [(j+Drows,j,1) for j in range(Prows)]

        self._store_matrix_list(newList, newShape)
    



def getCohomDimP(op1, op2, opP, rep_ind):
    """
    Computes the cohomology of the isotypical component corresponding to the projector P of the complex.
    :param op1: The operator corresponding to the differential D.
    : type op1: GraphOperator 
    :param op2: The operator corresponding to the differential DD, the cohomology is kerD/imDD.
    : type op2: GraphOperator
    :param opP: The projection operator corresponding to the isotypical component (i.e. im P).
    : type op1: GraphOperator  
    : return : The 
    """
    # tt = CHairyGraphComplex.ContractEdgesGO.generate_operator(
    #     n_vertices, n_loops, n_hairs, even_edges)
    # tu = CHairyGraphComplex.ContractEdgesGO.generate_operator(
    #     n_vertices+1, n_loops, n_hairs, even_edges)
    # symmp1 = CHairyGraphComplex.SymmProjector.generate_operator(
    #     n_vertices, n_loops, n_hairs, even_edges, rep_ind)

    D1 = op1.get_matrix()
    D2 = op2.get_matrix()
    # C = D2*D1
    opP.build_matrix(ignore_existing_files=True)
    P1 = opP.get_matrix()
    print("matrices loaded")

    D1P = D1*P1
    D2P = P1*D2
    print("computing ranks....")

    # diff = D1*D2
    # diffs = sum(abs(c) for cc in diff.columns() for c in cc)
    # print(diffs)
    
    isocomp_dim = P1.rank()
    r1 = D1P.rank()
    r2 = D2P.rank()
    # print(isocomp_dim, r1, r2)
    cohomdim = isocomp_dim - r1-r2
    n = opP.domain.get_n()
    part = Partitions(n)[rep_ind]
    rep_dim = symmetrica.charvalue(part, [1 for j in range(n)])
    # if cohomdim > 0:
    #     print("Cohomology found:  ", "even_edges" if even_edges else "odd_edges", ", h=", n_hairs, ", l=", n_loops, ", vertices=", n_vertices,
    #           " (degree ", n_vertices+n_loops -
    #           1, "), partition=", part,  ", invpartition=", part.conjugate(),
    #           ", multiplicity=", cohomdim/rep_dim, ", cohomdim=", cohomdim)
    # return isocomp_dim - r1-r2
    return cohomdim / rep_dim


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

