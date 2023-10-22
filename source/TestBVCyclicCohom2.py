import BVCyclic
import OrdinaryGraphComplex

import GenerateLatexTables

def cohom_formatted_forested_top(D1, D2, Dc2, use_Dc2_rank=None, iso_dict=None, use_D1domain = True):
    vs = D1.get_domain() if use_D1domain else D2.get_target()
    if not vs.is_valid():
        return "-"
    if not vs.exists_basis_file():
        return "?v"
    d = vs.get_dimension()

    r1 = 0
    r2 = 0
    rc2 = 0

    if Dc2.is_valid():
        if use_Dc2_rank is not None:
            if use_Dc2_rank == "?":
                return "?d"
            else:
                rc2 = int(use_Dc2_rank)
        elif Dc2.exists_rank_file():
            rc2 = Dc2.get_matrix_rank()
        else:
            return "?c2"

    if D1.is_valid():
        if D1.exists_rank_file():
            r1 = D1.get_matrix_rank()
        else:
            return "?1"
    if D2.is_valid():
        if D2.exists_rank_file():
            r2 = D2.get_matrix_rank()
        else:
            return "?2"

    # exact or not?
    r_str = "" 

    # iso string
    cohomdim = d+rc2-r1-r2
    iso_str = ""

    return str(cohomdim) + r_str + iso_str


def create_ordinarycyclic_cohom_table(v_range, l_range):
    s = ""

    header =  [str(v) for v in v_range]
    print(header)
    for even_edges in [False]:
        data = []
        for l in l_range:
            print(
                str(l),
                [ cohom_formatted_forested_top(
                    BVCyclic.ContractReconnectBiOM.generate_operator(v, l),
                    BVCyclic.ContractReconnectBiOM.generate_operator(v+1, l),
                    BVCyclic.ReconnectEdgesGO.generate_operator(v-2, l),
                    use_D1domain=False
                )  for v in v_range])
    return s

create_ordinarycyclic_cohom_table(range(4, 22), range(3, 12))


op1 = BVCyclic.ContractReconnectBiOM.generate_operator(10,6)
op2 = BVCyclic.ContractReconnectBiOM.generate_operator(11,6)
vs = op2.get_target()
opc = BVCyclic.ReconnectEdgesGO.generate_operator(8,6)
ovs1 = BVCyclic.GOneVS(9,6)
ovs2 = BVCyclic.GOneVS(8,6)

print(vs.get_dimension(), op1.get_matrix_rank(), op2.get_matrix_rank(), opc.get_matrix_rank())
print(ovs1.get_dimension(), ovs2.get_dimension())