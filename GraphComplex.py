from abc import ABCMeta, abstractmethod


class GraphComplex(object):
    __metaclass__ = ABCMeta

    def __init__(self, op_collection):
        self.op_collection = op_collection

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def get_cohomology_plot_path(self):
        pass

    @abstractmethod
    def plot_cohomology_dim(self):
        pass

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.op_collection.vs_collection.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                       progress_bar=progress_bar)

    def build_matrix(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.op_collection.build_matrix(ignore_existing_files=ignore_existing_files,
                                        n_jobs=n_jobs, progress_bar=progress_bar)

    def compute_rank(self, ignore_existing_files=True, n_jobs=1):
        self.op_collection.compute_rank(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs)

