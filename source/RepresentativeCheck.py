import GraphOperator
import GraphVectorSpace
from sage.all import *
from abc import ABCMeta, abstractmethod

class DummyVSOneMore(GraphVectorSpace.VectorSpace):
    def __init__(self, vs):
        self.vs = vs
        super(DummyVSOneMore, self).__init__()

    def get_dimension(self):
        return self.vs.get_dimension() + 1

    def is_valid(self):
        return self.vs.is_valid()

class RepresentativeCheck(GraphOperator.OperatorMatrix):
    """Class and method for checking whether a cocycle represents a non-zero cohomology class.
    """

    __metaclass__ = ABCMeta


    def __init__(self, op1, op2, name : str, primaryvs : GraphVectorSpace.GraphVectorSpace = None):
        self.op1 : GraphOperator.OperatorMatrix = op1
        self.op2 : GraphOperator.OperatorMatrix = op2
        self.name = name
        self.primaryvs = primaryvs
        # we only support that the graphs live in a single GraphVectorspace
        # by default, this is computed from op1
        if not primaryvs:
            if isinstance(self.op1.domain, GraphVectorSpace.SumVectorSpace):
                self.primaryvs = self.op1.domain.get_vs_list()[0]
            else:
                self.primaryvs = self.op1.domain
        print(str(self.primaryvs))
        super(RepresentativeCheck, self).__init__(DummyVSOneMore(op2.domain), op2.target)

    def __str__(self):
        return self.name

    def is_valid(self):
        return self.op2.is_valid()

    @staticmethod
    def is_match(domain, target):
        return True

    @abstractmethod
    def generate_vector(self):
        """returns the vector of the cocycle as a list of pairs (G, coeff)"""
        pass

    def get_vector_g6(self):
        vlst = self.generate_vector()
        tmp = ( (self.primaryvs.graph_to_canon_g6(G), c) for (G,c) in vlst )
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
        
        return v

    def get_vector1(self):
        """Returns the vector in the domain of op1"""
        return self._get_vector_rel(self.op1.domain)
    
    def get_vector2(self):
        """Returns the vector in the target of op2"""
        return self._get_vector_rel(self.op2.target)

    def is_cocycle(self):
        if not self.op1.is_valid():
            print(self.name, ": is a cocycle (trivially).")
            return True
        v = vector(self.get_vector1())
        A = self.op1.get_matrix()
        w = A*v 
        norm = w.norm(p=1)
        if norm == 0:
            print(self.name, ": is a cocycle.")
            return True
        else :
            print(self.name, f": is NOT a cocycle (norm={norm}).")
            return True

    def build_matrix(self, ignore_existing_files=False, skip_if_no_basis=True, progress_bar=True, **kwargs):
        """concatenates the vector to the matrix of op2"""
        if not self.op2.is_valid():
            return
        if (not ignore_existing_files) and self.exists_matrix_file():
            return

        vlst = self.get_vector2()
        (mlist, (m,n)) = self.op2._load_matrix_list()
        # print("*******", m,n,mlist)
        newshape = (m+1,n)
        newmlist = mlist + [(m,k,val) for (k,val) in enumerate(vlst) if val != 0 ]

        self._store_matrix_list(newmlist, newshape)

    def is_vector_zero(self):
        v = self.get_vector1()
        return sum(abs(c) for c in v) == 0

    def is_cocycle_exact(self):
        if self.is_vector_zero():
            print(self.name, f": is exact (since zero)")
            return True

        if not self.op2.is_valid():
            print(self.name, f": is not exact (trivially)")
            return False

        r1 = self.op2.get_matrix_rank()
        r2 = self.get_matrix_rank()
        if r1 > r2 or r2 > r1+1:
            raise RuntimeError(f"is_cocycle_exact error: r2 must be in [r1,r1+1] (r1={r1}, r2={r2})")
        
        if r1 == r2:
            print(self.name, f": is exact ({r1}={r2})")
        else: 
            print(self.name, f": is not exact ({r1}<{r2})")
        return r1 == r2

    def checkit(self, **kwargs):
        self.build_matrix()
        self.compute_rank(**kwargs)
        self.is_cocycle()
        self.is_cocycle_exact()