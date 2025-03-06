import OrdinaryGraphComplex
import CHairyGraphComplex


# V = OrdinaryGraphComplex.OrdinaryGVS(8,7, False)
V = CHairyGraphComplex.CHairyGraphVS(0,0,2,False)
# V = CHairyGraphComplex.CHairyGraphVS(5,2,3,False)
# Compute basis
V.build_basis()
# Iterate over basis, displaying every vector (as g6 ascii)
for s in V.get_basis_g6():
    print(s)
# Plot the basis elements to images, open html page (temp/temp.html) with the images
V.display_basis_plots()