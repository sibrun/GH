import WOHairyGraphComplex2

# v,l,h,w = 0,1,0,1
# v,l,h,w = 2,2,1,2
# v,l,h,w,c = 2,2,1,1,1
# v,l,h,w,c = 2,2,2,2,2
# v,l,h,w,c = 1,2,0,1,1

# v,l,h,w,c = 2,1,2,2,3
# v,l,h,w = 3,5,3,10
# v,l,h,w,c = 2,2,2,2,3
# v,l,h,w,c = 0,0,0,0,2

maxc = 15
maxg = 10
maxn = 12
maxv = 3
maxw = 12

max_excess = 5

for c in range(1,maxc+1):
    for g in range(maxg+1):
        for n in range(maxn+1):
            if 3*g+2*n - 25 <= max_excess:
                for v in range(maxv+1):
                    for w in range(maxw+1):
                        if 3*g+2*n-2*w <= max_excess:
                            V = WOHairyGraphComplex2.WOHairyGraphPreVS(v,g,n,w,c)
                            if V.is_valid(): 
                                print(V)
                                if not V.exists_basis_file():
                                    if c==1:
                                        for CV in V.get_chairy_prerequisites():
                                            CV.build_basis()
                                    V.build_basis(ignore_existing_files=False)
                                    # V.build_basis(ignore_existing_files=True)

