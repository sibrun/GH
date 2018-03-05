from abc import ABCMeta, abstractmethod


class GraphComplex(object):
    __metaclass__ = ABCMeta

    def __init__(self, vector_space, grading, differential):
        self.vector_space = vector_space
        self.grading = grading
        self.differential = differential

    @abstractmethod
    def get_cohomology_plot_path(self):
        pass

    @abstractmethod
    def plot_cohomology_dim(self):
        pass

    def get_vector_space(self):
        return self.vector_space

    def get_grading(self):
        return self.grading

    def get_differential(self):
        return self.differential

    def build_basis(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.vector_space.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                      progress_bar=progress_bar)

    def build_grading(self):
        self.grading.build_grading()

    def build_matrix(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        self.differential.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                       progress_bar=progress_bar)

    def compute_rank(self, ignore_existing_files=True, n_jobs=1):
        self.differential.compute_rank(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs)
