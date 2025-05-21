import KneisslerGC
import GraphVectorSpace
import GraphOperator
import OrdinaryGraphComplex


def print_cohom(nloops, even_edges):
    # print the cohomology of the graph complex
    V = KneisslerGC.KneisslerGVS(nloops, 0, even_edges)
    op_full = KneisslerGC.ContractEdgesGO.generate_operator(nloops, 2, even_edges)
    op_compl = KneisslerGC.ContractEdgesGO.generate_operator(nloops, 3, even_edges)
    r_full = op_full.get_matrix_rank()
    r_compl = op_compl.get_matrix_rank()
    n = V.get_dimension()

    op_ord = OrdinaryGraphComplex.ContractEdgesGO.generate_operator(2*nloops-2, nloops, even_edges)
    if op_ord.exists_rank_file():
        r_ord = op_ord.get_matrix_rank()
        n_ord = op_ord.domain.get_dimension()
        dim_ord = n_ord - r_ord
    else:
        dim_ord = "?"

    print("nloops=", nloops, "\tee=", even_edges, ": Actual cohomdim: ", dim_ord, " estimated cohomdim: ", n - r_full+ r_compl, "(full: ", r_full, " compl:", r_compl, ")")


for i in range(5, 11):
    print_cohom(i, True)

print("")
for i in range(5, 11):
    print_cohom(i, False)