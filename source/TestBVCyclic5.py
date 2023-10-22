import BVCyclic


ovs1 = BVCyclic.GOneVS(9,6)
ovs2 = BVCyclic.GOneVS(8,6)

l = list(ovs1.get_generating_graphs2())
print(len(l))
print(ovs1.graph_to_canon_g6( l[25]))
ovs1.build_basis(ignore_existing_files=True)
print(ovs1.get_dimension())

# print(vs.get_dimension(), op1.get_matrix_rank(), op2.get_matrix_rank(), opc.get_matrix_rank())
# print(ovs1.get_dimension(), ovs2.get_dimension())