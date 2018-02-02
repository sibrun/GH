import os
from joblib import Parallel, delayed
import multiprocessing
import Profiling
import logging
import StoreLoad as SL
import numpy as np
import Display

#reload(SL)
#reload(Profiling)

if __name__=='__main__':
    xrange = range(5,10)
    yrange = range(4,10)
    vdict = {(5,4):1,(7,8):2,(8,9):3}

    Display.save_2_indices_plot(vdict,'x',xrange,'y',yrange,'t','./test.png')

'''
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
'''





