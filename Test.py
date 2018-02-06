import HairyGraphComplex as HGC
import NautyInterface as NI
reload(NI)
reload(HGC)
v=5
l=4
h=4
even_edges = True
even_hairs = True

hgc = HGC.HairyGVS(v, l, h, even_edges, even_hairs)

generatingList = hgc._generating_graphs()

for g in generatingList:
    g.plot().show()

print(len(generatingList))

