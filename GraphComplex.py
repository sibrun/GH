from abc import ABCMeta, abstractmethod
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

    def __init__(self, vs_list, op_list):
        self.vs_list = vs_list
        self.op_list = op_list
        self.info_file_path = self._set_info_file_path()
        self.cohomology = dict()

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def _set_info_file_path(self):
        pass

    def members_to_string(self):
        vector_space = ["%s: %s" % (str(vs),vs.get_info()) for vs in self.vs_list]
        operator = ["%s: %s" % (str(op),op.get_info()) for op in self.op_list]
        return (vector_space, operator)

    def store_member_info(self):
        (vector_space, operator) = self.members_to_string()
        cohomology = self.get_cohomology_info()
        LHL = [("----- Graph Complex -----", [str(self)]),("----- Vector Space -----", vector_space),("----- Operator -----", operator),("----- Cohomology -----", cohomology)]
        SH.store_list_of_header_lists(LHL, self.info_file_path)

    def build_basis(self, skip_existing_files=True):
        self.vs_list.sort(key=operator.methodcaller('get_work_estimate'))
        for vs in self.vs_list:
            if not skip_existing_files:
                vs.delete_matrix_file()
            vs.build_basis()
        self.vs_list.sort(key=operator.methodcaller('get_dimension'))

    def build_operator_matrix(self, skip_existing_files=True):
        self.op_list.sort(key=operator.methodcaller('get_work_estimate'))
        for op in self.op_list:
            if not skip_existing_files:
                op.delete_matrix_file()
            op.build_matrix()
        self.op_list.sort(key=operator.methodcaller('get_matrix_entries'))

    def square_zero_test(self, eps):
        succ = []  # holds pairs for which test was successful
        fail = []  # failed pairs
        triv = []  # pairs for which test trivially succeeded because at least one operator is the empty matrix
        inc = []  # pairs for which operator matrices are missing
        for op1 in self.op_list:
            tvs1 = op1.target
            for op2 in self.op_list:
                dvs2 = op2.domain
                if tvs1 == dvs2:
                    # A composable pair is found
                    p = (op1, op2)
                    if not (op1.valid and op2.valid):
                        triv.append(p)
                        continue
                    try:
                        M1 = op1.get_matrix()
                        M2 = op2.get_matrix()
                    except SH.NotBuiltError:
                        logging.warn("Cannot test square zero: Operator matrix not built for %s or %s" % (str(op1), str(op2)))
                        inc.append(p)
                        continue
                    if op1.is_trivial() or op2.is_trivial():
                        triv.append(p)
                        continue
                    if sparse.linalg.norm(M2 * M1) < eps:
                        succ.append(p)
                    else:
                        fail.append(p)
        (triv_l, succ_l, inc_l, fail_l) = (len(triv), len(succ), len(inc), len(fail))
        logging.info("Square zero test for %s: trivial success: %d, success: %d, inconclusive %d, failed %d pairs" % (str(self), triv_l, succ_l, inc_l, fail_l ))
        if inc_l:
            logging.warn("Square zero test for %s: inconclusive: %d paris" % (str(self), inc_l))
        for (op1, op2) in fail:
            logging.error("Square zero test for %s: failed for the pair %s, %s" % (str(self), str(op1), str(op2)))
        return (triv_l, succ_l, inc_l, fail_l)

    #Computes the cohomology, i.e., ker(D)/im(DD)
    def compute_cohomology(self, only_dim=True):
        self.cohomology.clear()
        for opD in self.op_list:
            dvsD = opD.domain
            for opDD in self.op_list:
                tvsDD = opDD.target
                if tvsDD == dvsD:
                    (CH, conclusive) = GraphComplex._cohomology(opD, opDD, only_dim=only_dim)
                    if not conclusive:
                        self.cohomology.update({dvsD: "inconclusive"})
                        continue
                    self.cohomology.update({dvsD: CH})

    def get_cohomology_info(self):
        cohomologyList = []
        for vs in self.vs_list:
            dim = self.cohomology.get(vs)
            if dim is not None:
                line = "%s: dimension: %s" % (str(vs), str(dim))
            else:
                line = "%s: dimension: unknown" % str(vs)
            cohomologyList.append(line)
        return cohomologyList

    @staticmethod
    def _cohomology(opD, opDD, only_dim = True):
        if not opD.valid:
            D = 0
        else:
            try:
                D = opD.get_matrix()
            except SH.NotBuiltError:
                logging.warn("Cannot compute cohomology: Operator matrix not built for %s " % str(opD))
                return (None, False)
            if D.getnnz() == 0:
                D = 0
        if not opDD.valid:
                DD = 0
        else:
            try:
                DD = opD.get_matrix()
            except SH.NotBuiltError:
                logging.warn("Cannot compute cohomology: Operator matrix not built for %s " % str(opDD))
                return (None, False)
            if DD.getnnz() == 0:
                DD = 0

        if D is 0 and DD is 0:
            cohomology = opD.domain
            dim = opD.domain.get_dimension()
        else:                                       #TODO: Implement Cohomology
            if D is not 0 and DD is 0:
                A = D
            elif D is 0 and DD is not 0:
                A = DD
            else:
                pass #A = D * DD
        cohomology = "bla"
        dim = "not implemented"
        if only_dim:
            return (dim, True)
        return (cohomology, True)

