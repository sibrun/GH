from abc import ABCMeta, abstractmethod
import os
from sage.all import  *
import GraphVectorSpace as GVS

class GraphOperator():
    __metaclass__ = ABCMeta
    @abstractmethod
    def get_file_name(self):
        """Retrieve the file name (and path) of the file storing the matrix."""
        pass

    @abstractmethod
    def get_unique_file_name(self):
        """Retrieve a unique file name for the matrix.
           This filename is used when interchanging files with other computers."""
        pass

    @abstractmethod
    def get_source(self):
        """Returns the GraphVectorSpace{S} on which the operator acts."""
        pass

    @abstractmethod
    def et_target(self):
        """Returns the GraphVectorSpace{T} in which the operator takes values."""
        pass

    @abstractmethod
    def operate_on(self):
        """For G::S a graph in the domain, returns a list of pairs (GG, x), GG::T graph
        in the target, x a number,
        such that (operator)(G) = sum x GG."""
        pass

    @abstractmethod
    def get_work_estimate(self):
        """Provides a rough estimate of the amount of work needed to create the operator file.
          (In arbitrary units)"""
        pass

    def is_valid(self):
        vs = self.get_source()
        tvs = self.get_target()
        return vs.is_valid() and tvs.is_valid()

    def createOperatorFile(self):
        """
        Creates the matrix file that holds the operator.
        The corresponding list files for source and target
        must exist when calling this function.
        """
        outFile = self.get_file_name()
        vs = self.get_source()
        tvs = self.get_target()

        colorData = tvs.get_color_counts()

        if not os.path.isfile(domainListFile):
            print( "Cannot create operator file: First create list file $domainListFile")
            return

        if not os.path.isfile(targetListFile):
            print("Cannot create operator file: First create list file $targetListFile")
            return

        print( "Creating File $outFile ..." )

        sss = readAllLines(inListFile)
        lst = [from_string(S, s) for s in sss] # list of source graphs

        src_count = length(lst)
        ll = readAllLines(tgtListFile)
        tgt_count = length(ll)

        println( "List files read ($src_count, $tgt_count graphs)..." )
        #println("Test $(length(sss)) $(length(lst)) __ $(length(ll))")
       # mat = op.computeEdgeMarkingDifferential(ggg, tgtListFile, nMarks, evenEdges)

        count = 0 # counts current index in dummy file
        entries = [[] for G in lst] # will hold matrix

        if src_count == 0 || tgt_count == 0
            # create empty file and return
            open(outFile,"w") do f
            end
            println("Wrote empty file.")
            return
        end

        # lookup g6 -> index in target vector space
        lookup = Dict{String,Int}( s => j for (j,s) in enumerate(ll))

        open(outFile,"w") do f
          for (i,G) in enumerate(lst)
            the_image = operate_on(self, G)
            for (GG,prefactor) in the_image
                # canonize and look up
                GGcan, to_canonical_p = get_canon(GG, colorData)

                GGcang6 = to_string(GGcan)
                #println("$GGcang6 <- $(to_string(GG)): $to_canonical_p  ___ $(invPermutation( to_canonical_p ))")
                if haskey(lookup, GGcang6)
                  sgn = get_perm_sign(tvs, GGcan, to_canonical_p)
                  write(f, "$i $(lookup[GGcang6]) $(sgn * prefactor)\n" )
                end
            end
          end
          # write matrix size
          write(f, "$(length(lst)) $tgt_count 0\n")
          write(f, "$(length(lst)) $tgt_count 0\n")
        end

        println("done")

    def _read