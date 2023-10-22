import BVCyclic







ovs1 = BVCyclic.GOneVS(9,6)
ovs2 = BVCyclic.GOneVS(8,6)
print(list(ovs1.get_generating_graphs()))
l = list(ovs1.get_generating_graphs())
print(len(l))
ll = list(ovs1.get_generating_graphs2())
print(len(ll))
G = l[500]
autom_list = G.automorphism_group().gens()
print(G.graph6_string(), ovs1._has_odd_automorphisms(G,autom_list))
ovs1.build_basis(ignore_existing_files=True)
print(ovs1.get_dimension())

# print(vs.get_dimension(), op1.get_matrix_rank(), op2.get_matrix_rank(), opc.get_matrix_rank())
# print(ovs1.get_dimension(), ovs2.get_dimension())