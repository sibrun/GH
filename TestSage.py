from sage.all import *
import OrdinaryGraphComplex as OGC

reload(OGC)

nVertices=9
nLoops=8
print("----------------------------------------------------------------")
ogc=OGC.OrdinaryGraphVectorSpace(nVertices, nLoops, evenEdges=False)

#print(ogc.is_valid())
if ogc.valid():
    #print(ogc.get_file_name())
    #print(ogc.get_svg_dir())
    #print(ogc.get_work_estimate())

    #graphList=ogc._generating_graphs()
    #for g in graphList:
        #print(ogc.canonical_g6(g))
    #set_random_seed(1)
    #p=Permutations(range(1,nVertices+1)).random_element()
    #p=Permutation([2,1,3,4])
    ogc.create_basis()
    print(ogc.get_dimension())

