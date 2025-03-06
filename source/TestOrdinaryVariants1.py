import OrdinaryVariants


for l in range(9):
    for v in range(18):
        for even_edges in [True, False]:
            V = OrdinaryVariants.OrdinaryGVSFull(v,l,even_edges)
            V.build_basis()
            W = OrdinaryVariants.OrdinaryGVSBridgeless(v,l,even_edges)
            W.build_basis()
            U = OrdinaryVariants.OrdinaryGVSTriconnected(v,l,even_edges)
            U.build_basis()

