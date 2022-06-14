import HairyMerkulovComplex


  

ignore_ex = False

HGC = HairyMerkulovComplex.HairyMerkulovGC(range(15), range(7), range(1,3), True, True, differentials=['contract'])
HGC.build_basis(ignore_existing_files=ignore_ex)
HGC.build_matrix(ignore_existing_files=ignore_ex)
HGC.compute_rank(sage="integer",ignore_existing_files=ignore_ex)
print("Merkulov cohomology:")
HGC.print_cohom()
print("Reference cohomology:")
HGC.print_cohom_reference()