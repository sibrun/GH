from sage.all import *
from math import sqrt
from joblib import Parallel, delayed


def f (a,b):
    return a *b


x = Parallel(n_jobs=2, backend="threading")(delayed(f)(i,2) for i in range(10))
print(x)






