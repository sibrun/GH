import HairyGraphBiComplex as HGBC
import multiprocessing as mp

if __name__ == "__main__":

    '''ignore_ex = True

    gc = HGBC.HairyBiGC(range(6, 12), range(-5, -4) , True, False)

    gc.build_basis(n_jobs=1, info_tracker=True, ignore_existing_files=ignore_ex)

    gc.build_matrix(n_jobs=1, info_tracker=True, ignore_existing_files=ignore_ex)

    gc.square_zero_test()

    gc.compute_rank(info_tracker=True, ignore_existing_files=ignore_ex)

    gc.plot_cohomology_dim()'''

    def f(x):
        print(x)

    iter_arg = [1,2,3,4]
    pool = mp.Pool(3)
    [pool.apply_async(f, args=(x,)) for x in iter_arg]
    pool.close()
    pool.join()