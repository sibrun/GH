from OrdinaryGraphComplex import OrdinaryGraphVectorSpace
from sage.all import *

ogc=OrdinaryGraphVectorSpace(nVertices=4, nLoops=3)

print(ogc.is_valid())
if ogc.is_valid():
    print(ogc.get_file_name())
    print(ogc.get_svg_dir())
    print(ogc.get_work_estimate())

print(ogc.get_generating_graphs())

#gl=ogc.get_generating_graphs()