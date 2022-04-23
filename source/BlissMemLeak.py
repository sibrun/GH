
import os, psutil

from sage.all import *

process = psutil.Process(os.getpid())

oldmem = process.memory_info().rss
for i in range(1000000):
    G = graphs.RandomGNM(10,20)
    canonG = G.canonical_label(algorithm='bliss') # put bliss or sage
    # canonG = G.canonical_label(algorithm='sage')

    if i%1000 == 0:
        print(f"graph count {i}, mem usage (Delta) {process.memory_info().rss - oldmem}")
        oldmem = process.memory_info().rss

