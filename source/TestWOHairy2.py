import WOHairyGraphComplex2

# v,l,h,w = 0,1,0,1
# v,l,h,w = 2,2,1,2
# v,l,h,w,c = 2,2,1,1,1
# v,l,h,w,c = 2,2,2,2,2
# v,l,h,w,c = 1,2,0,1,1

v,l,h,w,c = 0,1,2,2,3
# v,l,h,w,c = 2,2,2,2,3
# v,l,h,w,c = 0,0,0,0,2

V = WOHairyGraphComplex2.WOHairyGraphPreVS(v,l,h,w,c)
# print (V.is_valid())
# for s in V.get_basis_g6():
#     print(s)
for CV in V.get_chairy_prerequisites():
    print(CV)
    CV.build_basis(ignore_existing_files=False)

for VV in V.get_self_prerequisites():
    print(VV)
    VV.build_basis(ignore_existing_files=True)

print("prereq done")
V.build_basis(ignore_existing_files=True)

V.display_basis_plots()