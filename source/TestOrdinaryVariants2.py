import OrdinaryVariants


for l in range(10):
    for v in range(18):
        for even_edges in [True, False]:
            D1 = OrdinaryVariants.ContractEdgesGOFull.generate_operator(v,l,even_edges)
            D1.build_matrix()
            D1.compute_rank(sage="integer")
            D2 = OrdinaryVariants.ContractEdgesGOBridgeless.generate_operator(v,l,even_edges)
            D2.build_matrix()
            D2.compute_rank(sage="integer")
            D3 = OrdinaryVariants.ContractEdgesGOTriconnected.generate_operator(v,l,even_edges)
            D3.build_matrix()
            D3.compute_rank(sage="integer")


