
import BVCyclic
import OrdinaryGraphComplex
import GraphVectorSpace
import GraphOperator


if __name__ == "__main__":
    nr_jobs = 10
    max_loops = 8

    print(f"Building GOneVS bases using {nr_jobs} jobs ...")
    sumvs = GraphVectorSpace.SumVectorSpace( [ BVCyclic.GOneVS(v, l) for v in range(0,16) for l in range(0,max_loops+1) ] )

    sumvs.build_basis(n_jobs=nr_jobs)

    print("Finished computing bases.")

    # op_list = [ BVCyclic.ReconnectEdgesGO.generate_operator(v, l) for v in range(0,20) for l in range(0,9) ]
    # op_list = [ BVCyclic.ContractReconnectBiOM.generate_operator(v, l) for v in range(0,20) for l in range(0,8) ]

    # allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)

    theD = BVCyclic.ContractReconnectTopD(range(16), range(max_loops+1))

    theD.build_matrix(n_jobs=nr_jobs)

    print("Finished computing matrices.")

    print("computing ranks")
    theD.compute_rank(linbox="rational", n_jobs=nr_jobs)
    print("Finished")