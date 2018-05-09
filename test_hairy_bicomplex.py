import HairyGraphBiComplex as HGBC
import HairyGraphComplex as HGC

if __name__ == "__main__":

    ignore_ex = False

    n_jobs = 1

    even_e = False
    even_h = False
    deg_range = range(3, 14)
    h_min_range = range(-11, -1)

    gc = HGBC.HairyBiGC(deg_range, h_min_range, even_e, even_h)

    gc.build_basis(info_tracker=True, ignore_existing_files=ignore_ex, n_jobs=n_jobs)

    gc.build_matrix(info_tracker=True, ignore_existing_files=ignore_ex, n_jobs=n_jobs)

    gc.square_zero_test()

    gc.compute_rank(info_tracker=True, ignore_existing_files=ignore_ex, n_jobs=n_jobs)

    gc.plot_cohomology_dim()