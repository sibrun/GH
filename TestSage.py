from sage.all import *
import operator
import random
import logging

G = graphs.WheelGraph(5)

for (j, ee) in enumerate(G.edges(labels=False)):
    a, b = ee
    G.set_edge_label(a, b, j)

print(G.edges())

G.merge_vertices([3,4])

print(G.edges())


print(type(G)==type(Graph()))
print(type(Graph()))
