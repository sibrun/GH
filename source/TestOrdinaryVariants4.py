import OrdinaryVariants
import OrdinaryGraphComplex

def print_dim_and_eulerchar_kv(v_range, l_range, k_range, even_edges):
    for k in k_range:
        print(k,"-vertex-connected")
        for l in l_range:
            ds = [OrdinaryVariants.OrdinaryGVSKVconnected(v, l,k, even_edges).get_dimension()
                    for v in v_range]
            eul = sum([(1 if j % 2 == 0 else -1) *
                        d for j, d in enumerate(ds)])
            print("Dimensions ",
                    l, even_edges, ":", ds, "Euler", eul)

def print_cohomology_dim_kv(v_range, l_range, k_range, even_edges):
    ret = {}
    for k in k_range:
        print(k,"-vertex-connected")
        for l in l_range:
            cohomdict = {}
            for v in v_range:
                D1 = OrdinaryVariants.ContractEdgesGOKV.generate_operator(
                    v, l, k, even_edges)
                D2 = OrdinaryVariants.ContractEdgesGOKV.generate_operator(
                    v+1, l, k, even_edges)
                try:
                    d = OrdinaryVariants.OrdinaryGVSKVconnected(v, l, k, even_edges).get_dimension()
                    r1 = D1.get_matrix_rank()
                    r2 = D2.get_matrix_rank()
                    cohomdict[v] = d-r1-r2
                except:
                    pass

            print("Cohomology Dimensions ",
                    l, even_edges, ":", cohomdict)
            ret[(k,l)] = cohomdict
    return ret

def print_dim_and_eulerchar_ke(v_range, l_range, k_range, even_edges):
    for k in k_range:
        print(k,"-edge-connected")
        for l in l_range:
            ds = [OrdinaryVariants.OrdinaryGVSKEconnected(v, l,k, even_edges).get_dimension()
                    for v in v_range]
            eul = sum([(1 if j % 2 == 0 else -1) *
                        d for j, d in enumerate(ds)])
            print("Dimensions ",
                    l, even_edges, ":", ds, "Euler", eul)

def print_cohomology_dim_ke(v_range, l_range, k_range, even_edges):
    ret = {}
    for k in k_range:
        print(k,"-edge-connected")
        for l in l_range:
            cohomdict = {}
            for v in v_range:
                D1 = OrdinaryVariants.ContractEdgesGOKE.generate_operator(
                    v, l, k, even_edges)
                D2 = OrdinaryVariants.ContractEdgesGOKE.generate_operator(
                    v+1, l, k, even_edges)
                try:
                    d = OrdinaryVariants.OrdinaryGVSKEconnected(v, l, k, even_edges).get_dimension()
                    r1 = D1.get_matrix_rank()
                    r2 = D2.get_matrix_rank()
                    cohomdict[v] = d-r1-r2
                except:
                    pass

            print("Cohomology Dimensions ",
                    l, even_edges, ":", cohomdict)
            ret[(k,l)] = cohomdict
    return ret


maxl = 8
krange = list(range(1,5))

for even_edges in [True, False]:
    print("KV:")
    print_dim_and_eulerchar_kv(range(15), range(maxl+1), krange, even_edges)
    r1=print_cohomology_dim_kv(range(15), range(maxl+1), krange, even_edges)
    print("KE:")
    print_dim_and_eulerchar_ke(range(15), range(maxl+1), krange, even_edges)
    r2=print_cohomology_dim_ke(range(15), range(maxl+1), krange, even_edges)
    # print("Ordinary (biconnected):")
    # print_dim_and_eulerchar_ordinary(range(20), range(maxl+2), even_edges)
    # print("Triconnected:")
    # print_dim_and_eulerchar_triconnected(range(20), range(maxl+2), even_edges)
    # r3=print_cohomology_dim_triconnected(range(15), range(maxl+1), even_edges)
    # print("Differences:")
    # print_differences(r1,r2,r3)
