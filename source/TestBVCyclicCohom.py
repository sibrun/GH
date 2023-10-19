import BVCyclic
import OrdinaryGraphComplex

def combirank(op1, op2):
    # computes rank of blockmatrix with both operators
    A= None
    B=None
    if op1.is_valid() and op1.domain.get_dimension() > 0 and op1.target.get_dimension() > 0:
        A = op1.get_matrix()
    if op2.is_valid() and op2.domain.get_dimension() > 0 and op2.target.get_dimension() > 0:
        B = op2.get_matrix()
    if A and B:
        return A.augment(B).rank()
    elif A:
        return A.rank()
    elif B:
        return B.rank()
    else:
        return 0
    

def cohom_dim(nv, nl):
    opc1 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(nv+1,nl, False)
    opc2 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(nv,nl, False)
    vs = opc1.target
    opr1 = BVCyclic.ReconnectEdgesGO.generate_operator(nv-1,nl)
    opr2 = BVCyclic.ReconnectEdgesGO.generate_operator(nv-2,nl)

    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    rc2 = 0

    if opr2.is_valid():
        if opr2.exists_rank_file():
            rc2 = opr2.get_matrix_rank()
        else:
            return "?"
    r1 = combirank(opc1, opr1)
    r2 = combirank(opc2, opr2)

    cohomdim = d+rc2-r1-r2

    return str(cohomdim) 

def cohom_dim2(nv, nl):
    opc1 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(nv+1,nl, False)
    opc2 = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(nv,nl, False)
    vs = opc1.target
    opr1 = OrdinaryGraphComplex.DeleteEdgesGO.generate_operator(nv,nl+1, False)
    opr2 = OrdinaryGraphComplex.DeleteEdgesGO.generate_operator(nv-1,nl+1, False)

    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?v"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    rc2 = 0

    if opr2.is_valid() and opr2.domain.get_dimension()>0 and opr2.target.get_dimension()>0:
        if opr2.exists_rank_file():
            rc2 = opr2.get_matrix_rank()
        else:
            return "?o"
    r1 = combirank(opc1, opr1)
    r2 = combirank(opc2, opr2)

    cohomdim = d+rc2-r1-r2

    return str(cohomdim) 


ooo = OrdinaryGraphComplex.DeleteEdgesGO.generate_operator(7,9,False)
ooo.build_matrix()
ooo.compute_rank(sage = "integer")
# print(ooo.is_valid(), ooo.get_matrix())

for nl in range(8):
    for nv in range(16):
        op = BVCyclic.ReconnectEdgesGO.generate_operator(nv,nl)
        # op.build_matrix(ignore_existing_files=True)
        op.build_matrix(ignore_existing_files=False, progress_bar=True)
        op.compute_rank(sage="integer", ignore_existing_files=False)
        # op.compute_rank(sage="integer", ignore_existing_files=True)

# Test commutativity
# for nv in range(9):
#     for nl in range(8):
#         opr1 = BVCyclic.ReconnectEdgesGO.generate_operator(nv-1,nl)
#         opr2 = BVCyclic.ReconnectEdgesGO.generate_operator(nv-2,nl)
#         opc = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(nv, nl, False)
#         if opr1.is_valid() and opr2.is_valid() and opc.is_valid() and opc.domain.get_dimension()>0 and opc.target.get_dimension()>0:
#             print("Testing ",nv, nl, opc.domain.get_dimension(), opc.target.get_dimension())
#             R1 = opr1.get_matrix()
#             R2 = opr2.get_matrix()
#             C = opc.get_matrix()
#             print(R1.dimensions(), R1.rank(), R2.dimensions(), R2.rank(), C.dimensions(), C.rank())
#             CR1 = C*R1
#             print(CR1.rank())
#             A = CR1.augment(R2)
#             print(A.rank())
#             print(R2.rank())



# compute cohomology:

for nl in range(8):
    print("nl=", nl)
    for nv in range(16):
        print(nv, "->", cohom_dim(nv, nl), cohom_dim2(nv, nl))