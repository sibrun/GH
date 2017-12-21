from abc import ABCMeta, abstractmethod
import os
import operator
import itertools
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO

reload(GVS)
reload(GO)


class GraphComplex():
    __metaclass__ = ABCMeta
    def __init__(self, delete_old=False):
        self.delete_old = delete_old
        self.vs_list = []
        self.op_list = []

    @abstractmethod
    def create_vs(self):
        pass

    @abstractmethod
    def create_op(self):
        pass

    def build_basis(self):
        self.create_vs()
        self.vs_list.sort(key=operator.attrgetter('work_estimate'))
        for vs in self.vs_list:
            if self.delete_old:
                vs.delete_file()
            vs.build_basis()

    def build_operator(self):
        self.create_op()
        self.op_list.sort(key=operator.attrgetter('work_estimate'))
        for op in self.op_list:
            if self.delete_old:
                op.delete_file()
            op.build_matrix()