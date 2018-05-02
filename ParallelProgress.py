import multiprocessing as mp
from multiprocessing.dummy import Pool as ThreadPool
import Queue
from tqdm import tqdm
import Parameters


def parallel_common_progress(func, iter_arg, common_arg, n_jobs=1, progress_bar=False, desc=None):
    if n_jobs == 1:
        if progress_bar:
            miniters = max(1, int(len(iter_arg) / Parameters.pbar_steps))
            return [func(x, common_arg) for x in tqdm(iter_arg, miniters=miniters, desc=desc)]
        else:
            print(desc)
            return [func(x, common_arg) for x in iter_arg]

    else:
        manager = mp.Manager()
        if type(common_arg) is dict:
            common_arg = manager.dict(common_arg)
        elif type(common_arg) is list:
            common_arg = manager.list(common_arg)

        if progress_bar:
            total = len(iter_arg)
            miniters = max(1, int(total/Parameters.pbar_steps))

            def update(res):
                pbar.update()
                return res

            with tqdm(total=total, miniters=miniters, desc=desc) as pbar:
                pool = mp.Pool(n_jobs)
                results = [pool.apply_async(func, args=(x, common_arg), callback=update) for x in iter_arg]
                result = [x.get() for x in results]

            return result

        else:
            print(desc)
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
            ptqdm = ParallelTqdm(iter_arg, results)
            queue = ptqdm.get_queue()
            for (idx, x) in enumerate(iter_arg):
                pbar_info = (True, True, idx, queue)
                pool.apply_async(func, args=(x,), kwds=dict({'pbar_info': pbar_info}, **kwargs),
                                 callback=results.append)
            pool.close()
            ptqdm.tqdm()
            pool.join()

        else:
            pbar_info = (False, False, None, None)
            [pool.apply_async(func, args=(x,), kwds=dict({'pbar_info': pbar_info}, **kwargs))
             for x in iter_arg]
            pool.close()
            pool.join()


def parallel(func, iter_arg, n_jobs=1, **kwargs):
    if n_jobs == 1:
        for x in iter_arg:
            func(x, **kwargs)

    else:
        pool = mp.Pool(n_jobs)
        [pool.apply_async(func, args=(x, ), kwds=kwargs) for x in iter_arg]
        pool.close()
        pool.join()


def parallel_thread(func, iter_arg, n_jobs=1, **kwargs):
    if n_jobs == 1:
        for x in iter_arg:
            func(x, **kwargs)

    else:
        pool = ThreadPool(n_jobs)
        [pool.apply_async(func, args=(x, ), kwds=kwargs) for x in iter_arg]
        pool.close()
        pool.join()


class ParallelTqdm(object):
    def __init__(self, iter_arg, results):
        self.n_bars = len(iter_arg)
        self.results = results
        manager = mp.Manager()
        self.queue = manager.Queue()
        self.pbars = [None] * self.n_bars
        self.position = 0

    def get_queue(self):
        return self.queue

    def start(self, idx, desc, total):
        self.pbars[idx] = tqdm(desc=desc, position=self.position, total=total)
        self.position += 1

    def update(self, idx, step):
        self.pbars[idx].update(step)

    def stop(self, idx):
        pbar = self.pbars[idx]
        pbar.close()

    def close(self):
        for pbar in self.pbars:
            pbar.close()

    def update_bars(self, message):
        (idx, mes, v, desc) = message
        if mes == 'start':
            self.start(idx, desc, v)
        if mes == 'step':
            self.update(idx, v)
        if mes == 'stop':
            self.stop(idx)

    def tqdm(self):
        while len(self.results) != self.n_bars:
            try:
                self.update_bars(self.queue.get(timeout=Parameters.pbar_timeout))
            except Queue.Empty:
                continue


def parallel_progress_messaging(func, iter_arg, arg, pbar_info=False, desc=None):
    if pbar_info is False:
        progress_bar = False
        print(desc)
    else:
        (progress_bar, message, idx, queue) = pbar_info
    if progress_bar:
        total = len(iter_arg)
        miniters = max(4, int(total / Parameters.pbar_steps))

        if message:
            queue.put((idx, 'start', total, desc))
            it_count = 0
        else:
            pbar = tqdm(total=total, desc=desc, miniters=miniters)

    for x in iter_arg:
        func(x, arg)
        if progress_bar:
            if message:
                it_count += 1
                if it_count % miniters == 0:
                    queue.put((idx, 'step', miniters, None))
            else:
                pbar.update()

    if progress_bar:
        if message:
            remaining = total - int(total / miniters) * miniters
            queue.put((idx, 'step', remaining, None))
            queue.put((idx, 'stop', None, None))
        else:
            pbar.close()
