import multiprocessing as mp
from tqdm import tqdm
import Parameters


def map_parallel(func, iter_arg, common_arg, n_jobs=1, progress_bar=False, desc=None):
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


def parallel_progress(func, iter_arg, n_jobs=1, progress_bar=False, *args, **kwargs):
    if n_jobs == 1:
        if progress_bar:
            for x in iter_arg:
                pbar = tqdm()
                func(x, 0, *args, **kwargs)
        else:
            for x in iter_arg:
                func(x, None, *args, **kwargs)

    else:
        if progress_bar:
            pool = mp.Pool(n_jobs)
            pbar_positions = range(len(iter_arg))
            [pool.apply_async(func, args=(x, pbar_pos, args, kwargs)) for pbar_pos, x in enumerate(iter_arg)]
            pool.close()
            pool.join()

        else:
            pool = mp.Pool(n_jobs)
            [pool.apply_async(func, args=(x, None, args, kwargs)) for x in iter_arg]
            pool.close()
            pool.join()

