from abc import ABCMeta, abstractmethod
import GraphVectorSpace as GVS
import GraphOperator as GO


class GraphComplex(object):
    __metaclass__ = ABCMeta

    def __init__(self, vector_space, operators):
        self.vector_space = vector_space
        self.operators = operators

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
        self.vector_space.build_basis(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs,
                                      progress_bar=progress_bar)

    def build_matrix(self, ignore_existing_files=True, n_jobs=1, progress_bar=False):
        for op in self.operators:
            op.build_matrix(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs, progress_bar=progress_bar)

    def compute_rank(self, ignore_existing_files=True, n_jobs=1):
        for op in self.operators:
            op.compute_rank(ignore_existing_files=ignore_existing_files, n_jobs=n_jobs)


class BiGraphComplex(GraphComplex):
    def __init__(self, vector_space, grading1, grading2, operators):

        super(BiGraphComplex, self).__init__(vector_space, operators)
        self.grading1 = grading1
        self.grading2 = grading2
        self.bi_grading = GVS.BiGrading(grading1, grading2)

    @abstractmethod
    def __str__(self):
        pass

    @abstractmethod
    def get_cohomology_plot_path(self):
        pass

    @abstractmethod
    def plot_cohomology_dim(self):
        pass

    def build_bi_grading(self):
        self.grading1.build_grading()
        self.grading2.build_grading()
        self.bi_grading.build_grading()
