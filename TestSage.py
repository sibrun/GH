from sage.all import *
import OrdinaryGraphComplex as OGC

reload(OGC)

nVertices=4
nLoops=3
print("----------------------------------------------------------------")
ogc=OGC.OrdinaryGraphVectorSpace(nVertices, nLoops, evenEdges=True)

#print(ogc.is_valid())
if ogc.is_valid():
    #print(ogc.get_file_name())
    #print(ogc.get_svg_dir())
    #print(ogc.get_work_estimate())

    graphList=ogc.get_generating_graphs()
    #set_random_seed(1)
    #p=Permutations(range(1,nVertices+1)).random_element()
    #p=Permutation([2,1,3,4])
    ogc.createListFile()
    print(ogc.getDimension())


