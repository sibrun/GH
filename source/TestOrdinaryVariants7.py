
import OrdinaryVariants
import GraphVectorSpace
import GraphOperator

"""Builds bases of all hairy vector spaces in the "computable" range.
"""


if __name__ == "__main__":
    nr_jobs = 6
    print(f"Building all computable variant matrices using {nr_jobs} jobs ...")
    vs_lista = []
    op_list = []
    maxl = 11

    for even_e in [True, False]:

        # vs_listf = vs_listf + [OrdinaryVariants.OrdinaryGVSFull(v, l, even_e)
        #                         for v in range(21) for l in range(maxl+1)]
        # vs_lista = vs_lista + [OrdinaryVariants.OrdinaryGVSBridgeless(v, l, even_e)
        #                         for v in range(21) for l in range(maxl+1)]
        vs_lista = vs_lista + [OrdinaryVariants.OrdinaryGVSTriconnected(v, l, even_e)
                                for v in range(17,21) for l in range(maxl, maxl+1)]

        # op_list = op_list + [OrdinaryVariants.ContractEdgesGOFull.generate_operator(v, l, even_e)
        #                         for v in range(21) for l in range(maxl+1)]
        # op_list = op_list + [OrdinaryVariants.ContractEdgesGOBridgeless.generate_operator(v, l, even_e)
        #                         for v in range(21) for l in range(maxl+1)]
        op_list = op_list + [OrdinaryVariants.ContractEdgesGOTriconnected.generate_operator(v, l, even_e)
                                for v in range(18,21) for l in range(maxl, maxl+1)]

    sumvs = GraphVectorSpace.SumVectorSpace(vs_lista)

    allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)


    # sumvsf.build_basis(n_jobs=nr_jobs, ignore_existing_files=True)
    print("Building other vector spaces.")
    sumvs.build_basis(n_jobs=nr_jobs)
    # sumvsa.build_basis(n_jobs=nr_jobs, ignore_existing_files=True)

    # sumvs.build_basis(n_jobs=nr_jobs)
    print("Building matrices.")
    allop.build_matrix(n_jobs=nr_jobs)
    # allop.build_matrix(n_jobs=nr_jobs, ignore_existing_files=True)

    print("Finished computing variant matrices.")

    print("computing ranks")
    allop.compute_rank(linbox="rational", n_jobs=nr_jobs)
    # allop.compute_rank(sage="integer", n_jobs=nr_jobs)
    # allop.compute_rank(sage="integer", n_jobs=nr_jobs, ignore_existing_files=True)
    print("Finished")
