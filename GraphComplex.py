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





    def build(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar)
        self.build_operator_matrix(ignore_existing_files=ignore_existing_files,
                                   n_jobs=n_jobs, progress_bar=progress_bar)









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
