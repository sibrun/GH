import HairyGraphComplex
import GraphVectorSpace

"""Builds bases of all hairy vector spaces in the "computable" range.

    computable means complexity <=max_complexity.
"""

max_complexity = 11

if __name__ == "__main__":
    nr_jobs = 1
    print(f"Building all computable hairy bases using {nr_jobs} jobs ...")
    vs_list = []

    for even_e in [True, False]:
        for even_h in [True, False]:
            # vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 1, even_e, even_h)
            #                      for v in range(18) for l in range(max_complexity-1)]
            # vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 2, even_e, even_h)
            #                      for v in range(20) for l in range(max_complexity-1)]
            # vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 3, even_e, even_h)
            #                      for v in range(20) for l in range(max_complexity-1)]
            # vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 4, even_e, even_h)
            #                      for v in range(20) for l in range(max_complexity-2)]
            # vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 5, even_e, even_h)
            #                      for v in range(20) for l in range(max_complexity-3)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 6, even_e, even_h)
                                 for v in range(20) for l in range(max_complexity-3)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 7, even_e, even_h)
                                 for v in range(20) for l in range(max_complexity-3)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v, l, 8, even_e, even_h)
                                 for v in range(20) for l in range(max_complexity-3)]

    sumvs = GraphVectorSpace.SumVectorSpace(vs_list)

    sumvs.build_basis(n_jobs=nr_jobs)

    print("Finished computing hairy bases.")
