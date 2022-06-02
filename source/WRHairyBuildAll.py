import WRHairyGraphComplex
import GraphVectorSpace

"""Builds bases of all wrhairy vector spaces in the "computable" range.
"""


if __name__ == "__main__":
    nr_jobs = 10
    print(f"Building all computable wrhairy bases using {nr_jobs} jobs ...")
    vs_list = []

    for w in [1, 2]:
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,0,w)
                            for v in range(18) for l in range(10)]
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,1,w)
                            for v in range(18) for l in range(9)]
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,2,w)
                            for v in range(18) for l in range(8)]
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,3,w)
                            for v in range(18) for l in range(7)]
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,4,w)
                            for v in range(18) for l in range(5)]
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,5,w)
                            for v in range(18) for l in range(4)] 
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,6,w)
                            for v in range(18) for l in range(3)] 
        vs_list = vs_list + [WRHairyGraphComplex.WRHairyGraphVS(v,l,7,w)
                            for v in range(18) for l in range(2)]         

    sumvs = GraphVectorSpace.SumVectorSpace(vs_list)

    sumvs.build_basis(n_jobs = nr_jobs)

    print("Finished computing hairy bases.")