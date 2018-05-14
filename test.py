from sage.all import *
from sage.combinat.shuffle import ShuffleProduct

import BiColoredHairyGraphComplex as BCHGC

#shuffles = ShuffleProduct(list(range(0,2)), list(range(2,4)))
#for shuffle in shuffles:
    #print(shuffle)


sum_vs = BCHGC.BiColoredHairyGraphSumVS(range(3,8),range(0,5),range(0,5),range(0,5),False,True,True)
#sum_vs.build_basis(info_tracker=True, n_jobs=4)
dif = BCHGC.ContractEdgesD(sum_vs)
#dif.build_matrix(info_tracker=True, n_jobs=4)
#dif.square_zero_test()
#dif.compute_rank(info_tracker=True, n_jobs=4)
dif.plot_cohomology_dim()

