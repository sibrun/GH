""" This file contains auxiliary methods and classes for graph complexes 
that are equipped with a symmetric group action.
The typical example is a graph complex with N numbered hairs, with the 
symmetric group acting by permutation of the hair labels. 
"""

import itertools
from sage.all import *
import GraphVectorSpace
import GraphOperator
import GraphComplex
import Shared
import OrdinaryGraphComplex
import StoreLoad
import Parameters

class SymmetricProjectionOperator(GraphOperator.GraphOperator):
    """This abstract class encodes the projection operator to an isotypical component of the symmetric group action
        by permuting numbered hairs.
        Warning: The matrix stores not the projector, but projector * n! / rep_dimension, to have integral matrices.

        The main method to be implemented by the user is 

        Deriving 
    """

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

    def __init__(self, domain, n, rep_index):
        """Initialize the domain and target vector space of the contract edges graph operator.

        : param domain: Domain vector space of the operator. This is also the target.
        : type domain: GraphVectorSpace
        : param n: The symmetric group acting should be S_n
        : type n: int
        : param rep_index: The index of the representation in the list produced by Partitions(h).
        : type rep_index: int
        """
        self.rep_index = rep_index
        self.n = n

        super(SymmetricProjectionOperator, self).__init__(domain, domain)

        # pre-fill in representation and character
        self.rep_partition = Partitions(n)[rep_index]
        self.norm_char_perm = [(symmetrica.charvalue(self.rep_partition, p.cycle_type(
            )), self.vertex_permutation_from_permutation(p)) for p in Permutations(n)]

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



