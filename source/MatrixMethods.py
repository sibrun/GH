""" Contains auxiliary methods for handling sparse matrices.
"""
from typing import List, Tuple
import StoreLoad
import os

def matrix_stats(matrix_file):
    """Gathers and prints general information on the operatormatrix.
    Currently only used for testing.
    """
    print("Loading...: ", matrix_file)
    (lst,(m,n)) = load_sms_file(matrix_file)
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


def get_submatrix(lst: List[Tuple[int,int,int]], keeprow: List[bool], keepcol: List[bool]):
    """Submatrix of a matrix.
    keeprow and keepcol are boolean masks of length m and n with the matrix of size (m,n)
    indicating which rows and columns the submatrix should contain. """
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

def transpose_lst(lst: List[Tuple[int,int,int]] ):
    """Transpose of a matrix."""
    matrix_list = [(j, i, v) for (i, j, v) in lst]
    matrix_list.sort()
    return matrix_list

def _removerstep(lst: List[Tuple[int,int,int]], m,n, rankbias):
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

def load_sms_file(fname: str) -> Tuple[List[Tuple[int,int,int]] , Tuple[int, int]]:
    """Loads a matric from an sms file.
    Returns a pair of a matrix (list) and the matrix dimensions.
    """
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

def save_sms_file(lst, m, n, fname):
    """Saves a matrix of size m x n as file fname, in sms format."""
    (d, t) = m,n
    stringList = []
    stringList.append("%d %d %s" % (d, t, "M"))
    for (i, j, v) in lst:
        stringList.append("%d %d %d" % (i + 1, j + 1, v))
    stringList.append("0 0 0")
    StoreLoad.store_string_list(stringList,fname)


def precondition(mlst, mm, nn, ensure_m_greater_n = False):
    """Computes a preconditioned version of the operatormatrix, saving it to
    (original_matrix_filename)preconditioned_{rankbias}.txt
    rankbias is to be added to the rank to obtain the true rank, and is also returned.
    """
    # print("Loading...: ", str(op))
    # (lst,(m,n)) = op._load_matrix_list()
    m,n = (mm, nn)
    lst = mlst
    rankbias = 0
    for i in range(10):
        print(f"remover step {i}, rankbias {rankbias}...")
        lst, m, n, rankbias = _removerstep(lst, m, n, rankbias)

    # transpose to ensure we have fewer columns (-> linbox has better progress report)
    if ensure_m_greater_n and m< n:
        lst = transpose_lst(lst)
        m,n = (n,m)

    # save_sms_file(lst, m, n, op.get_matrix_file_path()+f".preconditioned_{rankbias}.txt")
    return (lst, (m,n), rankbias)

def precondition_file(matrix_file, ensure_m_greater_n = False):
    """preconditions, and saves matrix to (original_matrix_filename)preconditioned_{rankbias}.txt
    returns new filename and rankbias."""
    lst, (m,n) = load_sms_file(matrix_file)
    lst2, (m2,n2), rankbias = precondition(lst, m, n, ensure_m_greater_n=ensure_m_greater_n)
    precond_fname = matrix_file + f".preconditioned_{rankbias}.txt"
    save_sms_file(lst2, m2, n2, precond_fname)
    return (precond_fname, rankbias)
