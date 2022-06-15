import HairyMerkulovComplex

import HairyGraphComplex
  

ignore_ex = False
HGC = HairyGraphComplex.HairyGC(range(15), range(7), range(1,3), False, False, differentials=['contract'])
HGC.build_basis()
HGC.build_matrix()
HGC.compute_rank(sage="integer")

HGC = HairyMerkulovComplex.HairyMerkulovGC(range(15), range(7), range(1,3), False, False, differentials=['contract'])
HGC.build_basis(ignore_existing_files=ignore_ex)
HGC.build_matrix(ignore_existing_files=ignore_ex)
HGC.compute_rank(sage="integer",ignore_existing_files=ignore_ex)
print("Merkulov cohomology:")
HGC.print_cohom()
print("Reference cohomology:")
HGC.print_cohom_reference()