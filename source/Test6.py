import unittest
import itertools
import logging
import Log
import TestGraphComplex
import OrdinaryGraphComplex
from sage.all import *
import StoreLoad
import OrdinaryMerkulovComplex
import GraphOperator

# even_e = True
# OGC = OrdinaryMerkulovComplex.OrdinaryMerkulovGC(range(20), range(9), even_e, ['contract'])
# OGC.build_basis()
# OGC.build_matrix()
# OGC.compute_rank(sage="mod")

# OGC.print_cohom()


def precond_stats(op:GraphOperator.OperatorMatrix):
    print("Loading...: ", str(op))
    (lst,(m,n)) = op._load_matrix_list()
    print("gathering stats...")
    rcount = [0 for _ in range(m)]
    ccount = [0 for _ in range(n)]
    multis = 0

    previj = (-1,-1)
    for (i,j, v) in lst:
        rcount[i] += 1
        ccount[j] += 1
        if previj == (i,j):
            multis += 1
    
    print("evaluating...")
    zerorows = sum( 1 if j==0 else 0 for j in rcount)
    zerocols = sum( 1 if j==0 else 0 for j in ccount)
    onerows = sum( 1 if j==1 else 0 for j in rcount)
    onecols = sum( 1 if j==1 else 0 for j in ccount)
    
    print(f"Matrix:    {m} x {n}")
    print(f"Zero rows: {zerorows}")
    print(f"Zero cols: {zerocols}")
    print(f"One rows:  {onerows}")
    print(f"One cols:  {onecols}")
    print(f"multis:    {multis}")


precond_stats(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(12,10, True))