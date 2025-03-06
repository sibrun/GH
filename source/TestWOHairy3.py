import WOHairyGraphComplex2

# v,l,h,w = 0,1,0,1
# v,l,h,w = 2,2,1,2
# v,l,h,w,c = 2,2,1,1,1
# v,l,h,w,c = 2,2,2,2,2
# v,l,h,w,c = 1,2,0,1,1

# v,l,h,w,c = 2,1,2,2,3
v,l,h,w = 3,5,3,10
# v,l,h,w,c = 2,2,2,2,3
# v,l,h,w,c = 0,0,0,0,2

V = WOHairyGraphComplex2.WOHairyGraphVS(v,l,h,w)
# print (V.is_valid())
# for s in V.get_basis_g6():
#     print(s)

cprereq = [CV for VV in V.get_prerequisites() for CV in VV.get_chairy_prerequisites()]
cprereq = list(set(cprereq))

for CV in cprereq:
    print(CV)
    CV.build_basis(ignore_existing_files=False)

prereq = V.get_prerequisites()
prereq.sort(key=lambda x : x.n_comp)

prereq2 = prereq + [VVV for VV in prereq for VVV in VV.get_self_prerequisites() ]
prereq2 = list(set(prereq2))
prereq2.sort(key=lambda x : x.n_comp)

for VV in prereq2:
    print(VV)
    VV.build_basis(ignore_existing_files=True)

print("prereq done")
V.build_basis(ignore_existing_files=True)

V.display_basis_plots()