import OrdinaryVariants

for l in range(8,11):
    print(l, " loops:")
    V = OrdinaryVariants.OrdinaryGVSFull(8, 6, False)
    print 
    for G in V.get_basis():
        print(G.graph6_string(),  G.edge_connectivity())