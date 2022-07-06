import HairyGraphComplex
import GraphVectorSpace
import GraphOperator

"""Builds bases of all hairy vector spaces in the "computable" range.
"""


if __name__ == "__main__":
    nr_jobs = 10
    print(f"Building all computable hairy matrices using {nr_jobs} jobs ...")
    vs_list = []
    op_list = []

    for even_e in [True, False]:
        for even_h in [True, False]:

            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 1, even_e, even_h)
                                 for v in range(18) for l in range(10)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 2, even_e, even_h)
                                 for v in range(18) for l in range(9)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 3, even_e, even_h)
                                 for v in range(18) for l in range(7)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 4, even_e, even_h)
                                 for v in range(18) for l in range(6)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 5, even_e, even_h)
                                 for v in range(18) for l in range(6)]

            op_list = op_list + [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, 1, even_e, even_h)
                                 for v in range(18) for l in range(10)]
            op_list = op_list + [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, 2, even_e, even_h)
                                 for v in range(18) for l in range(9)]
            op_list = op_list + [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, 3, even_e, even_h)
                                 for v in range(18) for l in range(7)]
            op_list = op_list + [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, 4, even_e, even_h)
                                 for v in range(18) for l in range(6)]
            op_list = op_list + [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, 5, even_e, even_h)
                                 for v in range(18) for l in range(6)]

    sumvs = GraphVectorSpace.SumVectorSpace(vs_list)

    allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)

    # sumvs.build_basis(n_jobs=nr_jobs)
    allop.build_matrix(n_jobs=nr_jobs)

    print("Finished computing hairy matrices.")
