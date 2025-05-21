import KneisslerGC
import GraphVectorSpace
import GraphOperator
import Parameters


# remove if not needed!!!
# Parameters.canonical_label_algorithm = "bliss"


if __name__ == "__main__":
    nr_jobs = 15
    print(f"Building all computable variant matrices using {nr_jobs} jobs ...")
    vs_listb = []
    vs_lista = []
    op_list = []
    maxl = 10
    typesa = [0,1,2]
    typesb = [3]
    types_op = [2,3]

    for even_e in [True, False]:

        vs_lista = vs_lista + [KneisslerGC.KneisslerGVS(l, t, even_e)
                                for t in typesa for l in range(4,maxl+1)]
        vs_listb = vs_listb + [KneisslerGC.KneisslerGVS(l, t, even_e)
                                for t in typesb for l in range(4,maxl+1)]

        op_list = op_list + [KneisslerGC.ContractEdgesGO.generate_operator(l, t, even_e)
                                for t in types_op for l in range(5,maxl+1)]

    sumvsb = GraphVectorSpace.SumVectorSpace(vs_listb)
    sumvsa = GraphVectorSpace.SumVectorSpace(vs_lista)
    sumvs = GraphVectorSpace.SumVectorSpace(vs_lista + vs_listb)

    allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)

    print("Building a vector spaces.")
    # sumvsa.build_basis(n_jobs=nr_jobs)
    sumvsa.build_basis(n_jobs=nr_jobs, ignore_existing_files=True)
    print("Building b vector spaces.")
    # sumvsb.build_basis(n_jobs=nr_jobs)
    sumvsb.build_basis(n_jobs=nr_jobs, ignore_existing_files=True)

    # sumvs.build_basis(n_jobs=nr_jobs)
    print("Building matrices.")
    # allop.build_matrix(n_jobs=nr_jobs)
    allop.build_matrix(n_jobs=nr_jobs, ignore_existing_files=True)

    print("Finished computing variant matrices.")

    print("computing ranks")
    # allop.compute_rank(linbox="mod", n_jobs=nr_jobs)
    allop.compute_rank(sage="integer", n_jobs=nr_jobs, ignore_existing_files=True)
    # allop.compute_rank(linbox="mod", n_jobs=nr_jobs, ignore_existing_files=True)
    # allop.compute_rank(sage="integer", n_jobs=nr_jobs)
    # allop.compute_rank(sage="integer", n_jobs=nr_jobs, ignore_existing_files=True)
    print("Finished")
