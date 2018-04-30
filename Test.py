import HairyGraphBiComplex as HGBC
import multiprocessing as mp

if __name__ == "__main__":

    gc = HGBC.HairyBiGC(range(6, 12), -4, True, False)

    '''gc.build_basis(n_jobs=1, info_tracker=True, ignore_existing_files=False)

    gc.build_matrix(n_jobs=1, info_tracker=True, ignore_existing_files=False)

    gc.square_zero_test()

    gc.compute_rank(info_tracker=True, ignore_existing_files=False)'''

    gc.plot_cohomology_dim()