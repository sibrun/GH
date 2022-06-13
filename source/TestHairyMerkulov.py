import HairyMerkulovComplex

HGC = HairyMerkulovComplex.HairyMerkulovGC(range(15), range(5), range(1,3), True, True, differentials=['contract'])
HGC.build_basis()
HGC.build_matrix()
HGC.compute_rank(sage="integer")
HGC.print_cohom()