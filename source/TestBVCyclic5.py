
import BVCyclic
import OrdinaryGraphComplex
import GraphVectorSpace
import GraphOperator



if __name__ == "__main__":
    nr_jobs = 20
    max_loops = 9
    max_vert = 16
    ignore_ex = False

    sumvs = GraphVectorSpace.SumVectorSpace( [ OrdinaryGraphComplex.OrdinaryGVS(v,l,False) for v in range(max_vert+1) for l in range(max_loops+1) ] )

    op_list = [ BVCyclic.AddVReconnectEdgesGO.generate_operator(v,l) for v in range(max_vert+1) for l in range(max_loops+1) ]
    allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)
    allop.build_matrix(n_jobs=nr_jobs, ignore_existing_files=ignore_ex)


# ovs1 = BVCyclic.GOneVS(9,6)
# ovs2 = BVCyclic.GOneVS(8,6)
# print(list(ovs1.get_generating_graphs()))
# l = list(ovs1.get_generating_graphs())
# print(len(l))
# ll = list(ovs1.get_generating_graphs2())
# print(len(ll))
# G = l[500]
# autom_list = G.automorphism_group().gens()
# print(G.graph6_string(), ovs1._has_odd_automorphisms(G,autom_list))
# ovs1.build_basis(ignore_existing_files=True)
# print(ovs1.get_dimension())

# print(vs.get_dimension(), op1.get_matrix_rank(), op2.get_matrix_rank(), opc.get_matrix_rank())
# print(ovs1.get_dimension(), ovs2.get_dimension())
