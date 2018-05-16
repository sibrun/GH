from sage.all import *
from sage.combinat.shuffle import ShuffleProduct

import BiColoredHairyGraphComplex as BCHGC
import BiColoredHairyGraphBiComplex as bicomplex

#shuffles = ShuffleProduct(list(range(0,2)), list(range(2,4)))
#for shuffle in shuffles:
    #print(shuffle)

n_jobs = 4
ignore_ex = False

#gc = BCHGC.BiColoredHairyGC(range(5,9),range(0,5),range(0,6),range(0,6),False,True,True, ['contract', 'split'])

gc = bicomplex.BiColoredHaryContractSplitBiGC(range(4,13),range(-11,1),range(-10,1),False,True,True)

#gc.build_basis(info_tracker=True, n_jobs=n_jobs, ignore_existing_files=ignore_ex)
#gc.build_matrix(info_tracker=True, n_jobs=n_jobs, ignore_existing_files=ignore_ex)
#gc.square_zero_test()
#gc.test_pairwise_anti_commutativity()
#gc.compute_rank(info_tracker=True, n_jobs=n_jobs, n_primes=0, estimate=True, small_primes=True, ignore_existing_files=ignore_ex)
gc.plot_cohomology_dim(to_csv=True, x_plots=3)
