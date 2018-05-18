"""Provides parallel mapping of a function to an iterable argument.."""

import multiprocessing as mp


def parallel(func, iter_arg, n_jobs=1, **kwargs):
    """Maps the function func on the iterable iter_arg and executes it using n_jobs parallel processes.
    :param func: function object: Function to be maped on the ietrable argument.
    :param iter_arg: iterable: Iterable argument.
    :param n_jobs: positive int: Number of parallel processes.
    :param kwargs: Keyword arguments to be passed forward to the function.
    """
    if n_jobs == 1:
        for x in iter_arg:
            func(x, **kwargs)

    else:
        pool = mp.Pool(n_jobs)
        [pool.apply_async(func, args=(x, ), kwds=kwargs) for x in iter_arg]
        pool.close()
        pool.join()