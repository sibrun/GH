import OrdinaryVariants

l=6
for v in range(8,11):
    print(l, " loops, ",v," vertices:")
    V = OrdinaryVariants.OrdinaryGVSFull(v, l, False)
    print 
    for G in V.get_basis():
        print(G.graph6_string(),  G.edge_connectivity())