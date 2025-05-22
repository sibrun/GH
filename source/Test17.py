import OrdinaryGraphComplex
from sage.all import *
from sage.matrix.berlekamp_massey import berlekamp_massey

D = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(10, 6, False)

theprime = 3323
p=theprime

At = D.get_matrix()
At.change_ring(GF(theprime))

A = At.transpose()

m = A.nrows()
n = A.ncols()
v = vector(GF(p), n)
for i in range(n):
    v[i] = 1

u = vector(GF(p), n)
for i in range(n):
    u[i] = 1

def wiedemann(A, Nlen, p):
    """
    Wiedemann algorithm for computing the kernel of a matrix A over GF(p).
    n is the number of rows of A.
    p is the prime field GF(p).
    """
    # Create a random vector in GF(p)^n
    m = A.nrows()
    n = A.ncols()
    v = vector(GF(p), n)
    for i in range(n):
        v[i] = 1

    u = vector(GF(p), n)
    for i in range(n):
        u[i] = 1
    
    # Compute the product of A and v
    res = vector(GF(p), Nlen)
    curv = v
    res[0] = u * curv
    At = A.transpose()
    for j in range(Nlen-1):
        curv = At * (A * curv)
        res[j+1] = u * curv

    return res
    

print(A.nrows(), A.ncols())
print(A*v)
print(At*(A*v))


seq = wiedemann(A, 30, theprime)

print("Wiedemann sequence: ", seq)

# print(At*A*v)

bmres = berlekamp_massey([GF(theprime)(s) for s in seq])
print("Berlekamp-Massey result: ", bmres, " degree: ", bmres.degree())
