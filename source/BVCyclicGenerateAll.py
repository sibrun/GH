
import BVCyclic
import OrdinaryGraphComplex
import GraphVectorSpace
import GraphOperator


if __name__ == "__main__":
    nr_jobs = 20
    max_loops = 9
    max_vert = 16
    ignore_ex = False

    print(f"Building GOneVS bases using {nr_jobs} jobs ...")
    sumvs = GraphVectorSpace.SumVectorSpace( [ BVCyclic.GOneVS(v, l) for v in range(0,max_vert+1) for l in range(0,max_loops+1) ] )

    sumvs.build_basis(n_jobs=nr_jobs, ignore_existing_files=ignore_ex)

    print("Finished computing bases.")

    # op_list = [ BVCyclic.ReconnectEdgesGO.generate_operator(v, l) for v in range(0,20) for l in range(0,9) ]
    # op_list = [ BVCyclic.ContractReconnectBiOM.generate_operator(v, l) for v in range(0,20) for l in range(0,8) ]

    # allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)

    op_list = [ BVCyclic.ReconnectEdgesGO.generate_operator(v,l) for v in range(0,max_vert+1) for l in range(0,max_loops+1) ]
    allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)
    allop.build_matrix(n_jobs=nr_jobs, ignore_existing_files=ignore_ex)
    allop.compute_rank(sage="integer", n_jobs=nr_jobs, ignore_existing_files=ignore_ex)

    theD = BVCyclic.ContractReconnectTopD(range(max_vert+1), range(max_loops+1))

    theD.build_matrix(n_jobs=nr_jobs, ignore_existing_files=ignore_ex)

    print("Finished computing matrices.")

    print("computing ranks")
    theD.compute_rank(sage="integer", n_jobs=nr_jobs, ignore_existing_files=ignore_ex)
    # theD.compute_rank(linbox="rational", n_jobs=nr_jobs)
    print("Finished")
