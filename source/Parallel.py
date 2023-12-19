"""Provide parallel mapping of a function to an iterable argument.."""

import multiprocessing as mp


def parallel(func, iter_arg, n_jobs=1, **kwargs):
    """Map the function func on the iterable iter_arg and executes it using n_jobs parallel processes.
    :param func: Function to be mapped on the iterable argument.
    :type func: function object
    :param iter_arg: Iterable argument.
    :type iter_arg: iterable
    :param n_jobs: Number of parallel processes.
    :type n_jobs: int
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
