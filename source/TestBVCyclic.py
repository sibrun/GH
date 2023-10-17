import BVCyclic


# Test basis generation
# for nv in range(8):
#     for nl in range(8):
#         BVCyclic.GOneVS(nv,nl).build_basis(ignore_existing_files=True)

# vs = BVCyclic.GOneVS(5,6)
# vs.display_basis_plots()

# print(vs.is_valid())

# vs = BVCyclic.GOneVS(4,3)
# vs.display_basis_plots()

# Test Operator generation
for nv in range(4):
    for nl in range(4):
        op = BVCyclic.ReconnectEdgesGO.generate_operator(nv,nl)
        op.build_matrix(ignore_existing_files=True)