import GraphOperator
import GraphVectorSpace
from sage.all import *
from abc import ABCMeta, abstractmethod

class RepresentativeCheck(GraphOperator.OperatorMatrix):
    """Class and method for checking whether a cocycle represents a non-zero cohomology class.
    """

    __metaclass__ = ABCMeta


    def __init__(self, op1, op2):
        self.op1 : GraphOperator.OperatorMatrix = op1
        self.op2 : GraphOperator.OperatorMatrix = op2

        super().__init__(op1.domain, op1.target)

    @abstractmethod
    def generate_vector(self):
        """returns the vector of the cocycle as a list of pairs (G, coeff)"""
        pass

    def get_vector_g6(self):
        vlst = self.generate_vector()
        tmp = ( (self.op1.domain.graph_to_canon_g6(G), c) for (G,c) in vlst )
        return [(s, c*sgn) for ((s,sgn),c) in tmp]

    def _get_vector_rel(self, vs):
        """Returns the vector in the domain of op1"""
        domdim = vs.get_dimension()
        v = [0 for _ in range(domdim)]
        basis_dict = vs.get_g6_coordinates_dict()
        vlst = self.get_vector_g6()
        for (s,c) in vlst:
            if s in basis_dict:
                j = basis_dict[s]
                v[j] = v[j] + c
        
        return vlst

    def get_vector1(self):
        """Returns the vector in the domain of op1"""
        return self._get_vector_rel(self.op1.domain)
    
    def get_vector2(self):
        """Returns the vector in the target of op2"""
        return self._get_vector_rel(self.op2.target)

    def is_cocycle(self):
        v = vector(self.get_vector1())
        A = self.op1.get_matrix()
        w = A*v 
        return w.norm(p=1) == 0

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        """concatenates the vector to the matrix"""
        vlst = self.get_vector2()
        (mlist, (m,n)) = self.op2.get_matrix_list()
        newshape = (m+1,n)
        newmlist = mlist + [(m+1,k,val) for (k,val) in enumerate(vlst) if val != 0 ]

        self._store_matrix_list(newmlist, newshape)

    def is_cocycle_exact(self):
        r1 = self.op2.get_rank()
        r2 = self.get_rank()
        if r1 > r2 or r2 > r1+1:
            raise RuntimeError(f"is_cocycle_exact error: r2 must be in [r1,r1+1] (r1={r1}, r2={r2})")
        return r1 == r2