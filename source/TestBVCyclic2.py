import BVCyclic
import OrdinaryGraphComplex


# Test basis generation
# for nv in range(12):
#     for nl in range(9):
#         # BVCyclic.GOneVS(nv,nl).build_basis(ignore_existing_files=True)
#         BVCyclic.GOneVS(nv,nl).build_basis(ignore_existing_files=False)

# vs = BVCyclic.GOneVS(5,6)
# vs.display_basis_plots()

# print(vs.is_valid())

# vs = BVCyclic.GOneVS(4,3)
# vs.display_basis_plots()



# Test Operator generation
for nv in range(12):
    for nl in range(9):
        op = BVCyclic.AddVReconnectEdgesGO.generate_operator(nv,nl)
        # op.build_matrix(ignore_existing_files=True, progress_bar=True)
        op.build_matrix(ignore_existing_files=False, progress_bar=True)
        op.compute_rank(sage="integer", ignore_existing_files=False)
        # op.compute_rank(sage="integer", ignore_existing_files=True)

        # op2=BVCyclic.ReconnectEdgesGO.generate_operator(nv, nl)
        # op2.domain.build_basis(ignore_existing_files=True)
        # op2.build_matrix(ignore_existing_files=True, progress_bar=True)
        # op2.compute_rank(sage="integer", ignore_existing_files=True)

# ooo = OrdinaryGraphComplex.DeleteEdgesGO.generate_operator(6, 8, False)
# print(ooo.is_valid(), ooo.domain.is_valid(), ooo.target.is_valid() )
# print(ooo.domain.get_dimension(), ooo.target.get_dimension())


# display and compare ranks
for nv in range(12):
    for nl in range(8):
        
        gvs = OrdinaryGraphComplex.OrdinaryGVS(nv+1,nl,False)
        if not gvs.is_valid() or gvs.get_dimension() == 0:
            continue

        print("Checking nv+1,nl=",nv+1,nl, " dimension ", gvs.get_dimension())
        # the target vector space has nv+1 vertices and nl loops
        opr = BVCyclic.ReconnectEdgesGO.generate_operator(nv,nl)
        print(opr.target)
        r1 = opr.get_matrix_rank()
        opd = BVCyclic.AddVReconnectEdgesGO.generate_operator(nv, nl+1)
        opn = OrdinaryGraphComplex.DeleteEdgesGO.generate_operator(nv+1, nl+1, False)
        print(opd.target)
        print(opn.target)
        r2 = opd.get_matrix_rank()
        r3 = opn.get_matrix_rank()
        print("Ranks: ",r1,r2, r3)

        if opr.is_valid() and opd.is_valid() and gvs.get_dimension() > 0:
            A = opr.get_matrix()
            B = opd.get_matrix()
            B2 = opn.get_matrix()
            print(B.dimensions(), B2.dimensions())
            BB = B.augment(B2)
            print(BB.rank())
            print(A.dimensions(), BB.dimensions())
            C = A.augment(BB)
            # C = A.stack(B)
            print(C.rank())
