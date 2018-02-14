import multiprocessing as mp
import Queue
from tqdm import tqdm
import Parameters


def parallel_common_progress(func, iter_arg, common_arg, n_jobs=1, progress_bar=False, desc=None):
    if n_jobs == 1:
        if progress_bar:
            miniters = int(len(iter_arg) / Parameters.pbar_steps)
            return [func(x, common_arg) for x in tqdm(iter_arg, miniters=miniters, desc=desc)]
        else:
            return [func(x, common_arg) for x in iter_arg]

    else:
        manager = mp.Manager()
        if type(common_arg) is dict:
            common_arg = manager.dict(common_arg)
        elif type(common_arg) is list:
            common_arg = manager.list(common_arg)

        if progress_bar:
            total = len(iter_arg)
            miniters = int(total/Parameters.pbar_steps)

            def update(res):
                pbar.update()
                return res

            with tqdm(total=total, miniters=miniters, desc=desc) as pbar:
                pool = mp.Pool(n_jobs)
                results = [pool.apply_async(func, args=(x, common_arg), callback=update) for x in iter_arg]
                result = [x.get() for x in results]

            return result

        else:
            pool = mp.Pool(n_jobs)
            results = [pool.apply_async(func, args=(x, common_arg)) for x in iter_arg]
            return [x.get() for x in results]


def parallel_individual_progress(func, iter_arg, n_jobs=1, progress_bar=False, **kwargs):
    if n_jobs == 1:
        pbar_info = (progress_bar, False, None, None)
        for x in iter_arg:
            func(x, pbar_info, **kwargs)

    else:
        pool = mp.Pool(n_jobs)
        if progress_bar:
            results = []
            iter_arg_l = len(iter_arg)
            pbars = [None] * iter_arg_l
            manager = mp.Manager()
            queue = manager.Queue()

            for (idx, x) in enumerate(iter_arg):
                pbar_info = (True, True, idx, queue)
                pool.apply_async(func, args=(x, pbar_info, kwargs), callback=results.append)
            pool.close()

            while len(results) != iter_arg_l:
                try:
                    update_pbars(pbars, queue.get(timeout=Parameters.timeout))
                except Queue.Empty:
                    continue
            pool.join()

        else:
            pbar_info = (False, False, None, None)
            [pool.apply_async(func, args=(x,), kwds=dict({'parallel_progress_bar_info': pbar_info}, **kwargs))
             for x in iter_arg]
            pool.close()
            pool.join()


def update_pbars(pbars, message):
    (idx, mes, v, desc) = message
    if mes == 'start':
        pbars[idx] = tqdm(desc=desc, position=idx, total=v)
    if mes == 'step':
        pbars[idx].update(v)
    if mes == 'stop':
        pbar = pbars[idx]
        pbar.n = pbar.total
        pbar.close()


def parallel(func, iter_arg, n_jobs=1, **kwargs):
    if n_jobs == 1:
        for x in iter_arg:
            func(x, **kwargs)

    else:
        pool = mp.Pool(n_jobs)
        [pool.apply_async(func, args=(x, kwargs)) for x in iter_arg]
        pool.close()
        pool.join()
