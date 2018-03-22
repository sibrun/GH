from abc import ABCMeta, abstractmethod
import Display


class GraphComplex(object):
    __metaclass__ = ABCMeta

    def __init__(self, vector_space, grading, differential):
        self.vector_space = vector_space
        self.grading = grading
        self.differential = differential

    @abstractmethod
    def get_param_labels_ranges_tuple(self):
        pass

    @abstractmethod
    def get_cohomology_plot_path(self):
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

    def square_zero_test(self):
        self.differential.square_zero_test()

    def compute_rank(self, exact=False, n_primes=1, estimate=True, ignore_existing_files=True, n_jobs=1):
        self.differential.compute_rank(exact=exact, n_primes=n_primes, estimate=estimate,
                                       ignore_existing_files=ignore_existing_files, n_jobs=n_jobs)

    def plot_cohomology_dim(self):
        dim_dict = self.differential.get_cohomology_dim()
        plot_path = self.get_cohomology_plot_path()
        param_labels_ranges = self.get_param_labels_ranges_tuple()
        Display.plot_array(dim_dict, param_labels_ranges, plot_path)
