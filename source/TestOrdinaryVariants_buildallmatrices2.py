import OrdinaryVariants
import GraphVectorSpace
import GraphOperator

"""Builds bases of all hairy vector spaces in the "computable" range.
"""


if __name__ == "__main__":
    nr_jobs = 15
    print(f"Building all computable variant matrices using {nr_jobs} jobs ...")
    vs_lista = []
    op_list = []
    maxl = 9
    maxk = 6

    for even_e in [True, False]:
        for k_conn in range(1,maxk+1):

            vs_lista = vs_lista + [OrdinaryVariants.OrdinaryGVSKVconnected(v, l, k_conn, even_e)
                                    for v in range(20) for l in range(maxl+1)]
            vs_lista = vs_lista + [OrdinaryVariants.OrdinaryGVSKEconnected(v, l, k_conn, even_e)
                                    for v in range(20) for l in range(maxl+1)]

            op_list = op_list + [OrdinaryVariants.ContractEdgesGOKV.generate_operator(v, l, k_conn, even_e)
                                    for v in range(20) for l in range(maxl+1)]
            op_list = op_list + [OrdinaryVariants.ContractEdgesGOKE.generate_operator(v, l, k_conn, even_e)
                                    for v in range(20) for l in range(maxl+1)]


    sumvs = GraphVectorSpace.SumVectorSpace(vs_lista)

    allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)


    print("Building other vector spaces.")
    sumvs.build_basis(n_jobs=nr_jobs)

    # sumvs.build_basis(n_jobs=nr_jobs)
    print("Building matrices.")
    allop.build_matrix(n_jobs=nr_jobs)

    print("Finished computing variant matrices.")

    print("computing ranks")
    allop.compute_rank(linbox="rational", n_jobs=nr_jobs)
    # allop.compute_rank(sage="integer", n_jobs=nr_jobs)
    print("Finished")
