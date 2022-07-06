import CHairyGraphComplex
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
        vs_list = vs_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 1, even_e)
                             for v in range(18) for l in range(8)]
        vs_list = vs_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 2, even_e)
                             for v in range(18) for l in range(8)]
        vs_list = vs_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 3, even_e)
                             for v in range(18) for l in range(8)]
        vs_list = vs_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 4, even_e)
                             for v in range(18) for l in range(7)]
        vs_list = vs_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 5, even_e)
                             for v in range(18) for l in range(6)]

        op_list = op_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 1, even_e)
                             for v in range(18) for l in range(8)]
        op_list = op_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 2, even_e)
                             for v in range(18) for l in range(8)]
        op_list = op_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 3, even_e)
                             for v in range(18) for l in range(8)]
        op_list = op_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 4, even_e)
                             for v in range(18) for l in range(7)]
        op_list = op_list + [CHairyGraphComplex.CHairyGraphVS(v, l, 5, even_e)
                             for v in range(18) for l in range(6)]

    sumvs = GraphVectorSpace.SumVectorSpace(vs_list)

    # sumvs.build_basis(n_jobs=nr_jobs)
    allop = GraphOperator.OperatorMatrixCollection(sumvs, op_list)

    # sumvs.build_basis(n_jobs=nr_jobs)
    allop.build_matrix(n_jobs=nr_jobs)

    print("Finished computing chairy matrices.")
