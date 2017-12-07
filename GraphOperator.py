from abc import ABCMeta, abstractmethod
import os
from sage.all import  *

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
