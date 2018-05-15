from sage.all import *
from sage.combinat.shuffle import ShuffleProduct

import BiColoredHairyGraphComplex as BCHGC
import BiColoredHairyGraphBiComplex as bicomplex

#shuffles = ShuffleProduct(list(range(0,2)), list(range(2,4)))
#for shuffle in shuffles:
    #print(shuffle)

n_jobs = 7

#gc = BCHGC.BiColoredHairyGC(range(5,8),range(0,5),range(0,5),range(0,5),False,True,True, ['contract', 'split'])

gc = bicomplex.BiColoredHaryContractSplitBiGC(range(4,11),range(-9,-1),range(-9,-1),False,True,True)

gc.build_basis(info_tracker=True, n_jobs=n_jobs)
gc.build_matrix(info_tracker=True, n_jobs=n_jobs)
gc.square_zero_test()
gc.test_pairwise_anti_commutativity()
gc.compute_rank(info_tracker=True, n_jobs=n_jobs)
gc.plot_cohomology_dim()
