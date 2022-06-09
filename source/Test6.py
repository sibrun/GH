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
import ForestedGraphComplex

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


def removerstep(lst, m,n, rankbias):
    print("gathering stats...")
    rcount = [0 for _ in range(m)]
    ccount = [0 for _ in range(n)]
    for (i,j, v) in lst:
        rcount[i] += 1
        ccount[j] += 1
    
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

    print("Simplifying matrix, deleting one cols")
    delrow = [False for _ in range(m)]
    delcol = [False for _ in range(n)]
    newrankbias = rankbias
    for (i,j, v) in lst:
        if ccount[j] == 1 and v != 0:
            if not delrow[i] and not delcol[j]:
                newrankbias += 1
            delrow[i] = True
            delcol[j] = True
    # remove zero cols
    for (j,v) in enumerate(ccount):
        if v ==0:
            delcol[j] = True
    for (i,v) in enumerate(rcount):
        if v ==0:
            delrow[i] = True


    newm = m - sum(1 for b in delrow if b)
    newn = n - sum(1 for b in delcol if b)
    
    newrowinds = [i for (i,b) in enumerate(delrow) if not b]
    newcolinds = [j for (j,b) in enumerate(delcol) if not b]
    
    # dict from old index to new
    rowdict = { iold : inew for (inew, iold) in enumerate(newrowinds) }
    coldict = { iold : inew for (inew, iold) in enumerate(newcolinds) }

    print("creating new matrix...")
    newlst = [ (rowdict[i],coldict[j],v) for (i,j,v) in lst 
                                if not delrow[i] and not delcol[j] ]

    return (newlst, newm, newn, newrankbias)

def precondition(op:GraphOperator.OperatorMatrix):
    print("Loading...: ", str(op))
    (lst,(m,n)) = op._load_matrix_list()
    rankbias = 0
    for i in range(10):
        print(f"remover step {i}, rankbias {rankbias}...")
        lst, m, n, rankbias = removerstep(lst, m, n, rankbias)

    (d, t) = m,n
    stringList = []
    stringList.append("%d %d %s" % (d, t, GraphOperator.data_type))
    for (i, j, v) in lst:
        stringList.append("%d %d %d" % (i + 1, j + 1, v))
    stringList.append("0 0 0")
    StoreLoad.store_string_list(stringList, op.get_matrix_file_path()+f".preconditioned_{rankbias}.txt")


    # save matrix


# precond_stats(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(12,10, True))
# precond_stats(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(13,10, True))
# precond_stats(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(11,10, False))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,9,0, True))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,10,0, True))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,8,0, True))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,10,0, False))

precondition( ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,9,0, False) )
# precondition( OrdinaryGraphComplex.ContractEdgesGO.generate_operator(12,10, True) )