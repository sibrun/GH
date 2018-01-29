import os
from joblib import Parallel, delayed
import multiprocessing
import Profiling
import logging
import StoreLoad as SL

#reload(SL)
#reload(Profiling)

log_path = os.path.join('log', 'test.log')
SL.generate_path(log_path)

logging.basicConfig(filename=log_path, level=logging.WARN)
logging.warn("###################################\n" + "----- Start test -----")

manager = multiprocessing.Manager()
res1 = manager.list()

res2 = list()

@Profiling.profile('log')
def f (a, b, res):
    res.append(a * b)
    return a * b


x = Parallel(n_jobs=2)(delayed(f)(i, 2, res1) for i in range(10))


print(x)
print(res1)

