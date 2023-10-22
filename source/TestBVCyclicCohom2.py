import BVCyclic
import OrdinaryGraphComplex

import GenerateLatexTables


def create_ordinarycyclic_cohom_table(v_range, l_range):
    s = ""

    header =  [str(v) for v in v_range]
    print(header)
    for even_edges in [False]:
        data = []
        for l in l_range:
            print(
                str(l),
                [ GenerateLatexTables.cohom_formatted_forested_top(
                    BVCyclic.ContractReconnectBiOM.generate_operator(v, l),
                    BVCyclic.ContractReconnectBiOM.generate_operator(v+1, l),
                    BVCyclic.ReconnectEdgesGO.generate_operator(v-2, l),
                    use_D1domain=False
                )  for v in v_range])
    return s

create_ordinarycyclic_cohom_table(range(4, 22), range(3, 12))

