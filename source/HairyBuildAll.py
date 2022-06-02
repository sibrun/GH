import HairyGraphComplex
import GraphVectorSpace

"""Builds bases of all hairy vector spaces in the "computable" range.
"""


if __name__ == "__main__":
    nr_jobs = 10
    print(f"Building all computable hairy bases using {nr_jobs} jobs ...")
    vs_list = []

    for even_e in [True, False]:
        for even_h in [True, False]:
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,1,even_e,even_h)
                                for v in range(18) for l in range(10)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,2,even_e,even_h)
                                for v in range(18) for l in range(9)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,3,even_e,even_h)
                                for v in range(18) for l in range(7)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,4,even_e,even_h)
                                for v in range(18) for l in range(6)]
            vs_list = vs_list + [HairyGraphComplex.HairyGraphVS(v,l,5,even_e,even_h)
                                for v in range(18) for l in range(6)]    

    sumvs = GraphVectorSpace.SumVectorSpace(vs_list)

    sumvs.build_basis(n_jobs = nr_jobs)

    print("Finished computing hairy bases.")