import WOHairyGraphComplex2
import cProfile

v,l,h,w,c = 3,2,11,12,10
# v,l,h,w,c = 1,0,11,9,9
# v,l,h,w,c = 1,3,9,11,8
# v,l,h,w,c = 2,2,2,2,3
# v,l,h,w,c = 0,0,0,0,2

V = WOHairyGraphComplex2.WOHairyGraphPreVS(v,l,h,w,c)
#cProfile.run('V.get_chairy_prerequisites()')
#print("*******")
cProfile.run('V.build_basis(ignore_existing_files=True)')
