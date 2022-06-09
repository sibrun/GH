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
import os

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
    tworows = sum( 1 if j==2 else 0 for j in rcount)
    twocols = sum( 1 if j==2 else 0 for j in ccount)
    
    print(f"Matrix:    {m} x {n}")
    print(f"Zero rows: {zerorows}")
    print(f"Zero cols: {zerocols}")
    print(f"One rows:  {onerows}")
    print(f"One cols:  {onecols}")
    print(f"Two rows:  {tworows}")
    print(f"Two cols:  {twocols}")
    print(f"multis:    {multis}")


def get_submatrix(lst, keeprow, keepcol):
    print("Submatrix...")
    newm = sum(1 for b in keeprow if b)
    newn = sum(1 for b in keepcol if b)
    
    newrowinds = [i for (i,b) in enumerate(keeprow) if b]
    newcolinds = [j for (j,b) in enumerate(keepcol) if b]
    
    # dict from old index to new
    rowdict = { iold : inew for (inew, iold) in enumerate(newrowinds) }
    coldict = { iold : inew for (inew, iold) in enumerate(newcolinds) }

    newlst = [ (rowdict[i],coldict[j],v) for (i,j,v) in lst 
                                if keeprow[i] and keepcol[j] ]

    return (newlst, newm, newn)

def graphcheck(lst, m, n):
    rcount = [0 for _ in range(m)]
    ccount = [0 for _ in range(n)]
    for (i,j, v) in lst:
        rcount[i] += 1
        ccount[j] += 1
    
    # consider only columns with 2 nonzero entries
    cols = [c == 2 for c in ccount]
    lst2, m2, n2 =  get_submatrix(lst, [True for _ in range(m)], cols)
    print("Generating graph")
    # make graph
    G = Graph(m2)
    edges = []
    es = [-1 for _ in range(n2)] 
    for (i,j,v) in lst2:
        if es[j] >0:
            edges.append( (i, es[j],j) )
        else:
            es[j] = i
    G.add_edges(edges)
    print("Graph generated, finding ccs")

    # ccgs = G.connected_components_subgraphs()
    # # ccsizes = [len(a) for a in ccs]
    # ccinfo = [ (GG.num_verts(), GG.num_edges()) for GG in ccgs ]
    # print(ccinfo[0:100])

    ccgs = G.connected_components_subgraphs()
    print("ccs done...")
    for cc in ccgs:
        # get submatrix
        verts = cc.vertices()
        if len(verts) < 100:
            break
        iset = set(verts)
        jset = set(cc.edge_labels())
        lst3, m3, n3 = get_submatrix(lst2, [(i in iset) for i in range(m2)], 
                [(j in jset) for j in range(n2)] )
        print(f"Submatrix {m3}x{n3}")
        M = matrix(QQ,m3,n3, {(i,j) : v for (i,j,v) in lst3 } )
        r = M.rank()
        if r == min(m3,n3):
            print( f"OK, full rank {r}")
        else:
            print( f"NO, rank deficient {r}")




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
    tworows = sum( 1 if j==2 else 0 for j in rcount)
    twocols = sum( 1 if j==2 else 0 for j in ccount)
    
    print(f"Matrix:    {m} x {n}")
    print(f"Zero rows: {zerorows}")
    print(f"Zero cols: {zerocols}")
    print(f"One rows:  {onerows}")
    print(f"One cols:  {onecols}")
    print(f"Two rows:  {tworows}")
    print(f"Two cols:  {twocols}")

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

def load_sms_file(fname):
    if not os.path.isfile(fname):
            raise StoreLoad.FileNotFoundError(
                "Cannot load matrix, No matrix file found for %s: " % fname)
    stringList = StoreLoad.load_string_list(fname)
    (d, t, data_type) = stringList.pop(0).split(" ")
    shape = (d, t) = (int(d), int(t))
    
    tail = map(int, stringList.pop().split(" "))
    if not list(tail) == [0, 0, 0]:
        raise ValueError("%s: End line missing or matrix not correctly read from file"
                            % fname)
    matrix_list = []
    for line in stringList:
        (i, j, v) = map(int, line.split(" "))
        if i < 1 or j < 1:
            raise ValueError("%s: Invalid matrix index: %d %d" %
                                (fname, i, j))
        if i > d or j > t:
            raise ValueError("%s: Invalid matrix index outside matrix size:"
                                " %d %d" % (fname, i, j))
        matrix_list.append((i - 1, j - 1, v))
    return (matrix_list, shape)

def precondition(op:GraphOperator.OperatorMatrix):
    print("Loading...: ", str(op))
    (lst,(m,n)) = op._load_matrix_list()
    rankbias = 0
    for i in range(10):
        print(f"remover step {i}, rankbias {rankbias}...")
        lst, m, n, rankbias = removerstep(lst, m, n, rankbias)

    (d, t) = m,n
    stringList = []
    stringList.append("%d %d %s" % (d, t, "M"))
    for (i, j, v) in lst:
        stringList.append("%d %d %d" % (i + 1, j + 1, v))
    stringList.append("0 0 0")
    # StoreLoad.store_string_list(stringList, op.get_matrix_file_path()+f".preconditioned_{rankbias}.txt")


    # save matrix


# precond_stats(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(12,10, True))
# precond_stats(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(13,10, True))
# precond_stats(OrdinaryGraphComplex.ContractEdgesGO.generate_operator(11,10, False))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,9,0, True))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,10,0, True))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,8,0, True))
# precond_stats(ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,10,0, False))

# precondition( ForestedGraphComplex.ContractUnmarkTopBiOM.generate_operator(7,9,0, False) )
# precondition( OrdinaryGraphComplex.ContractEdgesGO.generate_operator(12,10, True) )

print("Loading file")
fff = "gh_data/data/forestedbl/odd_edges/bi_D_contract_unmark_top_7_9_0.txt.preconditioned_978418.txt"

lst, (m,n) = load_sms_file(fff)
print("graphckeck...")
graphcheck(lst, m, n)