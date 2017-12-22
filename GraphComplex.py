from abc import ABCMeta, abstractmethod
import os
import operator
import itertools
import logging
import scipy.sparse as sparse
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO
import Shared as SH

reload(GVS)
reload(GO)
reload(SH)


class GraphComplex():
    __metaclass__ = ABCMeta
    def __init__(self, skip_existing_files=False):
        self.skip_existing_files = skip_existing_files
        self.vs_list = []
        self.op_list = []
        self.file_path = self._set_file_path()

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def _set_file_path(self):
        pass

    @abstractmethod
    def create_vs(self):
        pass

    @abstractmethod
    def create_op(self):
        pass

    def members_to_string(self):
        vector_space = ["%s: %s" % (str(vs),vs.get_info()) for vs in self.vs_list]
        operator = ["%s: %s" % (str(op),op.get_info()) for op in self.op_list]
        return (vector_space, operator)

    def store_member_list(self):
        self.load_info_from_file()
        (vector_space, operator) = self.members_to_string()
        LHL = [("----- Vector Space -----", vector_space),("----- Operator -----", operator)]
        SH.store_list_of_header_lists(LHL, self.file_path)

    def build_basis(self):
        self.create_vs()
        self.vs_list.sort(key=operator.attrgetter('work_estimate'))
        for vs in self.vs_list:
            if self.skip_existing_files:
                vs.delete_file()
            vs.build_basis()

    def build_operator(self):
        self.create_op()
        self.op_list.sort(key=operator.attrgetter('work_estimate'))
        for op in self.op_list:
            if self.skip_existing_files:
                op.delete_file()
            op.build_matrix()

    def load_info_from_file(self):
        for vs in self.vs_list:
            vs.load_info_from_file()
        for op in self.op_list:
            op.load_info_from_file()

    def square_zero_test(self, eps):
        succ = []  # holds pairs for which test was successful
        fail = []  # failed pairs
        triv = []  # pairs for which test trivially succeeded because at least one operator is the empty matrix
        inc = []  # pairs for which operator matrices are missing
        for op1 in self.op_list:
            dvs1 = op1.domain
            for op2 in self.op_list:
                tvs2 = op2.target
                if dvs1 == tvs2:
                    # A composable pair is found
                    p = (op1, op2)
                    if not (op1.valid and op2.valid):
                        triv.append(p)
                        continue
                    try:
                        D1 = op1.get_matrix()
                        D2 = op2.get_matrix()
                    except SH.NotBuiltError:
                        logging.warn("Cannot test square zero: Operator matrix not built for% %s or %s" % (str(op1), str(op2)))
                        inc.append(p)
                        print('not built')
                        continue
                    if op1.is_trivial() or op2.is_trivial():
                        triv.append(p)
                        continue
                    if sparse.linalg.norm(D2 * D1) < eps:
                        succ.append(p)
                    else:
                        fail.append(p)
        (triv_l, succ_l, inc_l, fail_l) = (len(triv), len(succ), len(inc), len(fail))
        logging.info("Square zero test for %s: trivial success: %d, success: %d, inconclusive %d, failed %d pairs" % (str(self), triv_l, succ_l, inc_l, fail_l ))
        logging.warn("Square zero test for %s: inconclusive: %d, failed: %d pairs" % (str(self), inc_l, fail_l ))

        for (op1, op2) in fail:
            logging.error("Square zero test for %s: failed for the pair %s, %s" % (str(self), str(op1), str(op2)))

        return (triv_l, succ_l, inc_l, fail_l)
