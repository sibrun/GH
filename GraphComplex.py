from abc import ABCMeta, abstractmethod
import logging
import operator
import itertools
import pandas
import StoreLoad as SL
import Display
import ParallelProgress as PP


class GraphComplex():
    __metaclass__ = ABCMeta

    def __init__(self, vs_list, op_list):
        self.vs_list = vs_list
        self.op_list = op_list

    @abstractmethod
    def get_type(self):
        pass

    @abstractmethod
    def get_params_range(self):
        pass

    @abstractmethod
    def get_params_names(self):
        pass

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def get_info_file_path(self):
        pass

    @abstractmethod
    def get_cohomology_plot_path(self):
        pass

    @abstractmethod
    def get_cohomology_dim(self):
        pass

    @abstractmethod
    def plot_cohomology_dim(self):
        pass

    def plot_info(self):
        vsList = []
        opList = []
        for vs in self.vs_list:
            vsList.append(list(vs.get_params()) + list(vs.get_info()))
        for op in self.op_list:
                opList.append(list(op.get_params()) + list(op.get_info()))
        vsColumns = list(self.get_params_names()) + ['valid', 'dimension']
        opColumns = list(self.get_params_names()) + ['valid', 'shape', 'entries', 'rank']
        vsTable = pandas.DataFrame(data=vsList, columns=vsColumns)
        opTable = pandas.DataFrame(data=opList, columns=opColumns)
        vsTable.sort_values(by=['valid', 'dimension'], inplace=True, na_position='last')
        vsTable.reset_index()
        opTable.sort_values(by=['valid', 'entries'], inplace=True, na_position='last')
        opTable.reset_index()
        Display.display_pandas_dfs([vsTable, opTable], self.get_info_file_path())

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.plot_info()
        self.sort_vs()
        PP.parallel_individual_progress(self._build_single_basis, self.vs_list, n_jobs=n_jobs, progress_bar=progress_bar,
                                        ignore_existing_files=ignore_existing_files)

    def _build_single_basis(self, vs, pbar_info, ignore_existing_files=True):
        vs.build_basis(pbar_info, ignore_existing_files=ignore_existing_files)

    def build_operator_matrix(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.plot_info()
        self.sort_op()
        for op in self.op_list:
            op.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar)

    def build(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar)
        self.build_operator_matrix(ignore_existing_files=ignore_existing_files,
                                   n_jobs=n_jobs, progress_bar=progress_bar)

    def sort_vs(self, work_estimate=True):
        if work_estimate:
            self.vs_list.sort(key=operator.methodcaller('get_work_estimate'))
        else:
            self.vs_list.sort(key=operator.methodcaller('get_sort_value'))

    def sort_op(self, work_estimate=True):
        if work_estimate:
            self.op_list.sort(key=operator.methodcaller('get_work_estimate'))
        else:
            self.op_list.sort(key=operator.methodcaller('get_sort_value'))

    def square_zero_test(self, eps):
        succ = []  # holds pairs for which test was successful
        fail = []  # failed pairs
        triv = []  # pairs for which test trivially succeeded because at least one operator is the empty matrix
        inc = []  # pairs for which operator matrices are missing
        for (op1, op2) in itertools.product(self.op_list, self.op_list):
            if op2.matches(op1):
                # A composable pair is found
                p = (op1, op2)
                if not (op1.valid and op2.valid):
                    triv.append(p)
                    continue
                try:
                    M1 = op1.get_matrix()
                    M2 = op2.get_matrix()
                except SL.FileNotFoundError:
                    logging.warn("Cannot test square zero: "
                                 "Operator matrix not built for %s or %s" % (str(op1), str(op2)))
                    inc.append(p)
                    continue
                if op1.is_trivial() or op2.is_trivial():
                    triv.append(p)
                    continue
                if sum(map(abs, (M2 * M1).list())) < eps:
                    succ.append(p)
                else:
                    fail.append(p)
        (triv_l, succ_l, inc_l, fail_l) = (len(triv), len(succ), len(inc), len(fail))
        logging.warn("Square zero test for %s: trivial success: "
                     "%d, success: %d, inconclusive: %d, failed: %d pairs" % (str(self), triv_l, succ_l, inc_l, fail_l ))
        if inc_l:
            logging.warn("Square zero test for %s: inconclusive: %d paris" % (str(self), inc_l))
        for (op1, op2) in fail:
            logging.error("Square zero test for %s: failed for the pair %s, %s" % (str(self), str(op1), str(op2)))
        return (triv_l, succ_l, inc_l, fail_l)

    def compute_ranks(self, ignore_existing_files=True, n_jobs=1):
        self.plot_info()
        self.sort_op(work_estimate=False)
        PP.parallel(self._compute_single_rank, self.op_list, n_jobs=n_jobs, ignore_existing_files=ignore_existing_files)

    def _compute_single_rank(self, op, ignore_existing_files=True):
        op.compute_rank(ignore_existing_files=ignore_existing_files)

    #Computes the cohomology, i.e., ker(D)/im(DD)
    def get_general_cohomology_dim_dict(self):
        cohomology_dim = dict()
        for (opD, opDD) in itertools.product(self.op_list, self.op_list):
            if opD.matches(opDD):
                dim = opD.cohomology_dim(opDD)
                cohomology_dim.update({opD.domain: dim})
        return cohomology_dim

    def get_cohomology_info(self, cohomology_dim=None):
        if cohomology_dim is None:
            return []
        cohomologyList = []
        for vs in self.vs_list:
            dim = cohomology_dim.get(vs)
            if dim is not None:
                line = "%s: %s" % (str(vs), str(dim))
                cohomologyList.append(line)
        return cohomologyList
