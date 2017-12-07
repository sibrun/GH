from abc import ABCMeta, abstractmethod
import os
from sage.all import  *
import GraphVectorSpace as GVS

class GraphOperator():
    __metaclass__ = ABCMeta
    @abstractmethod
    def file_name(self):
        """Retrieve the file name (and path) of the file storing the matrix."""
        pass

    @abstractmethod
    def domain(self):
        """Returns the GraphVectorSpace on which the operator acts."""
        pass

    @abstractmethod
    def target(self):
        """Returns the GraphVectorSpace in which the operator takes values."""
        pass

    @abstractmethod
    def operate_on(self):
        """For G::S a graph in the domain, returns a list of pairs (GG, x), GG::T graph
        in the target, x a number,
        such that (operator)(G) = sum x GG."""
        pass

    @abstractmethod
    def work_estimate(self):
        """Provides a rough estimate of the amount of work needed to create the operator file.
          (In arbitrary units)"""
        pass

    def valid(self):
        return self.domain().valid() and self.target().valid()

    def create_operator_matrix(self):
        """
        Creates the matrix file that holds the operator.
        The corresponding list files for source and target
        must exist when calling this function.
        """
        fileName = self.file_name()
        domain = self.domain()
        target = self.target()

        try:
            domainBasis = domain.basis()
        except GVS.NotBuiltError:
            raise GVS.NotBuiltError("Cannot bild operator matrix: First build basis of the domain")
        try:
            targetBasis = target.basis()
        except GVS.NotBuiltError:
            raise GVS.NotBuiltError("Cannot bild operator matrix: First build basis of the target")

        domainDim = len(domainBasis)
        targetDim = len(targetBasis)

        if domainDim == 0 and targetDim == 0:
            # create empty file and return
            open(fileName,"w").close()
            return

        matrix = []

        # lookup g6 -> index in target vector space
        lookup = Dict{String,Int}( s => j for (j,s) in enumerate(ll))

        f = open(fileName,"w")
        for (i,G) in enumerate(domainBasis):
            image = self.operate_on(G)
            for GG, prefactor in image:
                # canonize and look up
                GGcanon = domain.canonical(GG, domain.color_counts())

                GGcanon6 = to_string(GGcan)
                #println("$GGcang6 <- $(to_string(GG)): $to_canonical_p  ___ $(invPermutation( to_canonical_p ))")
                if haskey(lookup, GGcang6)
                  sgn = get_perm_sign(tvs, GGcan, to_canonical_p)
                  write(f, "$i $(lookup[GGcang6]) $(sgn * prefactor)\n" )
          # write matrix size
          write(f, "$(length(lst)) $tgt_count 0\n")
          write(f, "$(length(lst)) $tgt_count 0\n")
        end

        println("done")

    def _read