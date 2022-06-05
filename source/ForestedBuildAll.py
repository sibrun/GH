import ForestedGraphComplex
import GraphVectorSpace

"""Builds bases of all hairy vector spaces in the "computable" range.
"""


if __name__ == "__main__":
    nr_jobs = 10
    print(f"Building all computable forested bases and operators using {nr_jobs} jobs ...")
    for even_e in [False, True]:
        FGC = ForestedGraphComplex.ContractUnmarkTopD(range(8), range(16), range(1), even_e)
        FGC.build_basis()
        FGC.build_matrix(n_jobs=nr_jobs)
        
        FGC = ForestedGraphComplex.ContractUnmarkTopD(range(7), range(16), range(1,2), even_e)
        FGC.build_basis()
        FGC.build_matrix(n_jobs=nr_jobs)

        FGC = ForestedGraphComplex.ContractUnmarkTopD(range(6), range(16), range(2,3), even_e)
        FGC.build_basis()
        FGC.build_matrix(n_jobs=nr_jobs)

        FGC = ForestedGraphComplex.ContractUnmarkTopD(range(5), range(16), range(3,4), even_e)
        FGC.build_basis()
        FGC.build_matrix(n_jobs=nr_jobs)

        FGC = ForestedGraphComplex.ContractUnmarkTopD(range(4), range(16), range(4,5), even_e)
        FGC.build_basis()
        FGC.build_matrix(n_jobs=nr_jobs)

        FGC = ForestedGraphComplex.ContractUnmarkTopD(range(3), range(16), range(5,6), even_e)
        FGC.build_basis()
        FGC.build_matrix(n_jobs=nr_jobs)

    # for even_e in [True, False]:
    #     for even_h in [True, False]:
    #         vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,1,even_e,even_h)
    #                             for v in range(18) for l in range(10)]
    #         vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,2,even_e,even_h)
    #                             for v in range(18) for l in range(9)]
    #         vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,3,even_e,even_h)
    #                             for v in range(18) for l in range(7)]
    #         vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,4,even_e,even_h)
    #                             for v in range(18) for l in range(6)]
    #         vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,5,even_e,even_h)
    #                             for v in range(18) for l in range(6)]    



        print("Finished computing forested bases and operators.")