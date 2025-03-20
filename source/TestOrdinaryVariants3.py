import OrdinaryVariants
import OrdinaryGraphComplex

def print_dim_and_eulerchar_full(v_range, l_range, even_edges):
    for l in l_range:
        ds = [OrdinaryVariants.OrdinaryGVSFull(v, l, even_edges).get_dimension()
                for v in v_range]
        eul = sum([(1 if j % 2 == 0 else -1) *
                    d for j, d in enumerate(ds)])
        print("Dimensions ",
                l, even_edges, ":", ds, "Euler", eul)

def print_cohomology_dim_full(v_range, l_range, even_edges):
    ret = {}
    for l in l_range:
        cohomdict = {}
        for v in v_range:
            D1 = OrdinaryVariants.ContractEdgesGOFull.generate_operator(
                v, l, even_edges)
            D2 = OrdinaryVariants.ContractEdgesGOFull.generate_operator(
                v+1, l, even_edges)
            try:
                d = OrdinaryVariants.OrdinaryGVSFull(v, l, even_edges).get_dimension()
                r1 = D1.get_matrix_rank()
                r2 = D2.get_matrix_rank()
                cohomdict[v] = d-r1-r2
            except:
                pass

        print("Cohomology Dimensions ",
                l, even_edges, ":", cohomdict)
        ret[l] = cohomdict
    return ret


def print_dim_and_eulerchar_bridgeless(v_range, l_range, even_edges):
    for l in l_range:
        ds = [OrdinaryVariants.OrdinaryGVSBridgeless(v, l, even_edges).get_dimension()
                for v in v_range]
        eul = sum([(1 if j % 2 == 0 else -1) *
                    d for j, d in enumerate(ds)])
        print("Dimensions ",
                l, even_edges, ":", ds, "Euler", eul)

def print_cohomology_dim_bridgeless(v_range, l_range, even_edges):
    ret = {}
    for l in l_range:
        cohomdict = {}
        for v in v_range:
            D1 = OrdinaryVariants.ContractEdgesGOBridgeless.generate_operator(
                v, l, even_edges)
            D2 = OrdinaryVariants.ContractEdgesGOBridgeless.generate_operator(
                v+1, l, even_edges)
            try:
                d = OrdinaryVariants.OrdinaryGVSBridgeless(v, l, even_edges).get_dimension()
                r1 = D1.get_matrix_rank()
                r2 = D2.get_matrix_rank()
                cohomdict[v] = d-r1-r2
            except:
                pass

        print("Cohomology Dimensions ",
                l, even_edges, ":", cohomdict)
        ret[l] = cohomdict
    return ret

def print_cohomology_dim_bridgeless_tex(v_range, l_range, even_edges):
    ret = {}
    print( " & "+ " & ".join(str(v) for v in v_range) , "\\\\ \\hline")
    for l in l_range:
        print(l, end=" ")
        for v in v_range:
            D1 = OrdinaryVariants.ContractEdgesGOBridgeless.generate_operator(
                v, l, even_edges)
            D2 = OrdinaryVariants.ContractEdgesGOBridgeless.generate_operator(
                v+1, l, even_edges)

            VV = OrdinaryVariants.OrdinaryGVSBridgeless(v, l, even_edges)
            d = VV.get_dimension()
            if VV.is_valid():
                try:
                    r1 = D1.get_matrix_rank() if D1.get_domain().is_valid() and D1.get_domain().get_dimension()>0 and D1.get_target().is_valid() and D1.get_target().get_dimension()>0 else 0
                    r2 = D2.get_matrix_rank() if D2.get_domain().get_dimension()>0 and D2.get_target().get_dimension()>0 else 0
                    cdim = str(d-r1-r2)
                except:
                    cdim="?"
            else:
                cdim = "-"
            print("& ", cdim, end=" ")
        print("\\\\")
        print("\\hline")

    return ret
   
def print_dim_and_eulerchar_triconnected(v_range, l_range, even_edges):
    for l in l_range:
        ds = [OrdinaryVariants.OrdinaryGVSTriconnected(v, l, even_edges).get_dimension()
                for v in v_range]
        eul = sum([(1 if j % 2 == 0 else -1) *
                    d for j, d in enumerate(ds)])
        print("Dimensions ",
                l, even_edges, ":", ds, "Euler", eul)

def print_cohomology_dim_triconnected(v_range, l_range, even_edges):
    ret = {}
    for l in l_range:
        cohomdict = {}
        for v in v_range:
            D1 = OrdinaryVariants.ContractEdgesGOTriconnected.generate_operator(
                v, l, even_edges)
            D2 = OrdinaryVariants.ContractEdgesGOTriconnected.generate_operator(
                v+1, l, even_edges)
            try:
                d = OrdinaryVariants.OrdinaryGVSTriconnected(v, l, even_edges).get_dimension()
                r1 = D1.get_matrix_rank()
                r2 = D2.get_matrix_rank()
                cohomdict[v] = d-r1-r2
            except:
                pass

        print("Cohomology Dimensions ",
                l, even_edges, ":", cohomdict)
        ret[l] = cohomdict
    return ret

def print_dim_and_eulerchar_ordinary(v_range, l_range, even_edges):
    for l in l_range:
        ds = [OrdinaryGraphComplex.OrdinaryGVS(v, l, even_edges).get_dimension()
                for v in v_range]
        eul = sum([(1 if j % 2 == 0 else -1) *
                    d for j, d in enumerate(ds)])
        print("Dimensions ",
                l, even_edges, ":", ds, "Euler", eul)

def print_differences(r1,r2,r3):
    for l in r1:
        d1 = r1[l]
        d2 = r2[l]
        d3 = r3[l]
        print("Differences for l=",l)
        for v in d1:
            if d1[v] != d2[v] or d1[v] != d3[v]:
                print(v,":", d1[v], d2[v], d3[v])


   

maxl = 10
maxv = 20

for even_edges in [True, False]:
    print("Full:")
    print_dim_and_eulerchar_full(range(maxv+1), range(maxl+1), even_edges)
    r1=print_cohomology_dim_full(range(maxv+1), range(maxl+1), even_edges)
    print("Bridgeless:")
    print_dim_and_eulerchar_bridgeless(range(maxv+1), range(maxl+1), even_edges)
    r2=print_cohomology_dim_bridgeless(range(maxv+1), range(maxl+1), even_edges)
    # print("Ordinary (biconnected):")
    # print_dim_and_eulerchar_ordinary(range(20), range(maxl+1), even_edges)
    print("Triconnected:")
    print_dim_and_eulerchar_triconnected(range(maxv+1), range(maxl+1), even_edges)
    r3=print_cohomology_dim_triconnected(range(maxv+1), range(maxl+1), even_edges)
    print("Differences:")
    print_differences(r1,r2,r3)

print_cohomology_dim_bridgeless_tex(range(4,21), range(3,11), True)
print_cohomology_dim_bridgeless_tex(range(4,21), range(3,11), False)