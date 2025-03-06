import BVCyclic
import OrdinaryGraphComplex


# Test basis generation
for nv in range(12):
    for nl in range(9):
        # BVCyclic.GOneVS(nv,nl).build_basis(ignore_existing_files=True)
        BVCyclic.GOneVS(nv,nl).build_basis(ignore_existing_files=False)

# vs = BVCyclic.GOneVS(5,6)
# vs.display_basis_plots()

# print(vs.is_valid())

# vs = BVCyclic.GOneVS(4,3)
# vs.display_basis_plots()



# Test Operator generation
for nv in range(12):
    for nl in range(9):
        op = BVCyclic.ReconnectEdgesGO.generate_operator(nv,nl)
        # op.build_matrix(ignore_existing_files=True)
        op.build_matrix(ignore_existing_files=False)
        op.compute_rank(sage="integer", ignore_existing_files=False)
        # op.compute_rank(sage="integer", ignore_existing_files=True)

# ooo = OrdinaryGraphComplex.DeleteEdgesGO.generate_operator(6, 8, False)
# print(ooo.is_valid(), ooo.domain.is_valid(), ooo.target.is_valid() )
# print(ooo.domain.get_dimension(), ooo.target.get_dimension())

gvs66 = OrdinaryGraphComplex.OrdinaryGVS(6,6, False)
# gvs66.display_basis_plots()

# display and compare ranks
for nv in range(12):
    for nl in range(9):

        gvs = OrdinaryGraphComplex.OrdinaryGVS(nv+1,nl,False)
        if not gvs.is_valid() or gvs.get_dimension() == 0:
            continue

        print("Checking nv+1,nl=",nv+1,nl, " dimension ", gvs.get_dimension())
        # the target vector space has nv+1 vertices and nl loops
        opr = BVCyclic.ReconnectEdgesGO.generate_operator(nv,nl)
        print(opr.target)
        r1 = opr.get_matrix_rank()
        opd = OrdinaryGraphComplex.DeleteEdgesGO.generate_operator(nv+1, nl+1, False)
        print(opd.target)
        r2 = opd.get_matrix_rank()
        print("Ranks: ",r1,r2)

        if opr.is_valid() and opd.is_valid() and gvs.get_dimension() > 0:
            A = opr.get_matrix()
            B = opd.get_matrix()
            print(A.dimensions(), B.dimensions())
            C = A.augment(B)
            # C = A.stack(B)
            print(C.rank())
