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

maxc = 4
maxg = 10
maxn = 12
maxv = 5
maxw = 12

max_excess = 5


for g in range(maxg+1):
    for n in range(maxn+1):
        if 3*g+2*n - 25 <= max_excess:
            for v in range(maxv+1):
                for w in range(maxw+1):
                    if 3*g+2*n-2*w <= max_excess:
                        V = WOHairyGraphComplex2.WOHairyGraphVS(v,g,n,w)
                        if V.is_valid():
                            print(V)
                            V.build_basis(ignore_existing_files=False)
                            # V.build_basis(ignore_existing_files=True)

