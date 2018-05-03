import HairyGraphBiComplex as HGBC

if __name__ == "__main__":

    ignore_ex = True

    n_jobs = 8

    even_e = False
    even_h = False
    deg_range = range(5, 12)
    h_min_range = range(-9, -1)

    gc = HGBC.HairyBiGC(deg_range, h_min_range, even_e, even_h)

    gc.build_basis(info_tracker=True, ignore_existing_files=ignore_ex, n_jobs=n_jobs)

    gc.build_matrix(info_tracker=True, ignore_existing_files=ignore_ex, n_jobs=n_jobs)

    gc.square_zero_test()

    gc.compute_rank(info_tracker=True, ignore_existing_files=ignore_ex, n_jobs=n_jobs)

    gc.plot_cohomology_dim()
