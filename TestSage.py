from sage.all import *
from joblib import Parallel, delayed
import multiprocessing
import Profiling

reload(Profiling)

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




