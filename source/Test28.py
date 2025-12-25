# test whether harmonic reps have small entries

from sage.all import *
import random
import OrdinaryGraphComplex

nverts = 10
nloops = 8

op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(nverts,nloops,False)
# op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(10,6,False)

D = op.get_matrix()

op2 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(nverts+1,nloops,False)
# op = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(10,6,False)

DD = op2.get_matrix()


# find number of automorphisms
aut_factors = [GG.automorphism_group().cardinality() for GG in op.get_domain().get_basis()]
# print(aut_factors)

# change base ring of D to finite field
# D = D.change_ring(GF(100003))


def find_rank_support(D):
    m = D.nrows()
    n = D.ncols()
    r = D.rank()
    print("Matrix size: ", D.nrows(), "x", D.ncols(), "; rank: ", r)

    active_set = set(random.sample(range(n), r))
    unused_set = set(range(n)) - active_set
    l_active_set = list(active_set)
    Dr = D[:, l_active_set]
    rr = Dr.rank()

    while rr < r or len(l_active_set) >r:
        print("Current rank: ", rr, "; need ", r)
        print("Removing columns...")
        # reduce rank if possible
        to_delete = len(l_active_set)-rr
        KK = Dr.right_kernel_matrix()
        print("Kernel size: ", KK.nrows(), " x ", KK.ncols())
        pcols = list(KK.pivots())
        if len(pcols) != to_delete:
            print("Error, have to delete ", to_delete, " but only ", len(pcols), " pivots in kernel")
            return
        # remove first to_delete columns that have a pivot in KK
        for pc in pcols:
            col_index = l_active_set[pc]
            active_set.remove(col_index)
            unused_set.add(col_index)

        if rr == r:
            break

        print("Adding columns")
        # add  random columns from unused set
        to_add = r-rr+50
        l_unused_set = list(unused_set)
        random_add = random.sample(l_unused_set, to_add)
        for ac in random_add:
            unused_set.remove(ac)
            active_set.add(ac)
        l_active_set = list(active_set)
        Dr = D[:, l_active_set]
        rr = Dr.rank()
    print("Final active set size: ", len(active_set))
    return l_active_set




m = D.nrows()
n = D.ncols()
r = D.rank()
print("Matrix size: ", D.nrows(), "x", D.ncols(), "; rank: ", r)

rs = find_rank_support(D)
print("Found rank support of size ", len(rs))
print(rs)

# high_aut_indices = [i for i in range(n) if aut_factors[i] > 1]
# hal = len(high_aut_indices)
# print(len(high_aut_indices), " : ", high_aut_indices)
# # check for column sets to delete
# # size 1
# for k in range(0-30,180):
#     for l in range(5):
#         indices = random.sample(range(n), r+k)
#         # indices = high_aut_indices + random.sample(range(n), k)
#         D2 = D[:, indices]
#         rr = D2.rank()
#         print("k= ", k, "Rank: ", rr)


# lll = DD.nonzero_positions()

# for (i,j) in lll:
#     DD[i,j] *= aut_factors[i]



# lll = D.nonzero_positions()

# for (i,j) in lll:
#     D[i,j] *= aut_factors[j]



# C = D.stack(DD.transpose())

# KK = C.right_kernel_matrix()

# print(KK)
# print(KK.ncols(), KK.nrows(), len(aut_factors))

# v = [KK[0,j] for j in range(len(aut_factors))]

# print(v)

# print(gcd(v))


# V = D.right_kernel()

# print(V)
# p = 1009  # 4-digit prime
# p = 100003  # 6-digit prime
# D_finite = D.change_ring(GF(p))
# V_finite = D_finite.right_kernel()
# print(V_finite)