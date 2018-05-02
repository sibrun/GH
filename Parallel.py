import multiprocessing as mp

def parallel(func, iter_arg, n_jobs=1, **kwargs):
    if n_jobs == 1:
        for x in iter_arg:
            func(x, **kwargs)

    else:
        pool = mp.Pool(n_jobs)
        [pool.apply_async(func, args=(x, ), kwds=kwargs) for x in iter_arg]
        pool.close()
        pool.join()