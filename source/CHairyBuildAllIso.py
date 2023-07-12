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
    CGCs = []
    for even_e in [True, False]:
        CGCs.append( CHairyGraphComplex.CHairyGC(range(0, 18), range(8), range(2, 3), even_e, ['contract_iso']) )
        CGCs.append( CHairyGraphComplex.CHairyGC(range(0, 18), range(7), range(3, 4), even_e, ['contract_iso']) )
        CGCs.append(CHairyGraphComplex.CHairyGC(range(0, 18), range(6), range(4, 5), even_e, ['contract_iso']) )
        CGCs.append(CHairyGraphComplex.CHairyGC(range(0, 18), range(5), range(5, 6), even_e, ['contract_iso']) )

        
    print("Computing matrices")
    for CGC in CGCs:
        CGC.build_matrix(progress_bar=False, info_tracker=False, n_jobs=nr_jobs)

    print("computing ranks")
    for CGC in CGCs:
        CGC.compute_rank(linbox="mod",progress_bar=False, info_tracker=False, n_jobs=nr_jobs)
    print("Finished")
