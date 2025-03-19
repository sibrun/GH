import unittest
import logging
import Log
import TestGraphComplex
import WHairyGraphComplex
from sage.all import Graph

log_file = "WHGC_Unittest.log"


def check_graphs_vs_basis(GVS, w):
    # Takes a list of graphs and checks whether they are found in the basis
    ba = GVS.get_basis_g6()
    for HH in w:
        H = Graph(HH) if type(HH) is str else HH
        g6, sgn = GVS.graph_to_canon_g6(H)
        autom_list = H.automorphism_group(partition=GVS.get_partition()).gens()
        if GVS._has_odd_automorphisms(H, autom_list):
            print(g6, " has odd automorphisms")
        else:
            if g6 not in ba:
                print(g6, " not found in basis")
            else:
                print(g6, " exists with index ", ba.index(g6))


def DSquareTestSingle(n_vertices, n_loops, n_hairs, n_ws, j_to_pick=-1):
    tt = WHairyGraphComplex.ContractEdgesGO.generate_operator(n_vertices, n_loops, n_hairs, n_ws)
    tu = WHairyGraphComplex.ContractEdgesGO.generate_operator(n_vertices-1, n_loops, n_hairs, n_ws)
    D1 = tt.get_matrix()
    D2 = tu.get_matrix()
    C = D2*D1

    if j_to_pick < 0:
        for i in range(C.nrows()):
            for j in range(C.ncols()):
                if C[i, j] != 0:
                    print(i, j, C[i, j])
                    j_to_pick = j
        if j_to_pick < 0:
            print("success, squares to zero")
            return
        else:
            print("Does not square to zero, checking index ", j_to_pick)
    # print(D1)
    # print(D2)
    # print(C)
    ba0 = tt.domain.get_basis_g6()
    ba1 = tu.domain.get_basis_g6()
    ba2 = tu.target.get_basis_g6()

    G = Graph(ba0[j_to_pick])
    w = tt.operate_on(G)

    # check whether graphs are in basis
    for H, x in w:
        g6, sgn = tu.domain.graph_to_canon_g6(H)
        autom_list = H.automorphism_group(partition=tu.domain.get_partition()).gens()
        if tu.domain._has_odd_automorphisms(H, autom_list):
            print(g6, " has odd automorphisms")
        else:
            if g6 not in ba1:
                print(g6, " not found in basis ", " v=",x)
            else:
                print(g6, " exists at index ", ba1.index(g6), " v=",x)

    # compute D^2
    ww = [(HH, x*xx) for H, x in w for HH, xx in tu.operate_on(H)]
    wwd = {}
    for H, x in ww:
        g6, sgn = tu.target.graph_to_canon_g6(H)
        if g6 in wwd:
            wwd[g6] += x
        else:
            wwd[g6] = x
    print(wwd)
    nonzeroflag = False
    for g6, x in wwd.items():
        if x != 0:
            print("Nonzero entry: ", g6, x)
            nonzeroflag = True
    if not nonzeroflag:
        print("all entries zero, i.e., success.")

# v_range = range(1, 8)
# l_range = range(0, 7)
# h_range = range(1, 8)
# edges_types = [True, False]
# hairs_types = [True, False]

# tt = WHairyGraphComplex.WHairyGraphVS(4,4,1,2)
# tt.build_basis()
# #tt.plot_all_graphs_to_file()
# tt.display_basis_plots()

#tt = WHairyGraphComplex.WHairyGraphVS(4,3,0,1)
#print(tt.is_valid())
#tt.build_basis(ignore_existing_files=True)

#tt.plot_all_graphs_to_file()
#tt.display_basis_plots()


# WGC = WHairyGraphComplex.WHairyGC(range(0,11), range(5,7), range(0,4), range(2,3) , ['contract'])
WGC = WHairyGraphComplex.WHairyGC(range(6), range(2), range(4), range(1, 2), ['contract'])
# WGC = WHairyGraphComplex.WHairyGC(range(0,8), range(0,6), range(1,3), range(2,3) , ['contract'])

WGC.build_basis(progress_bar=False, info_tracker=False, ignore_existing_files=True)
WGC.build_matrix(progress_bar=False, info_tracker=False, ignore_existing_files=True)

# WGC.build_basis(progress_bar=False, info_tracker=False, ignore_existing_files=False)
# WGC.build_matrix(progress_bar=False, info_tracker=False, ignore_existing_files=False)

WGC.square_zero_test()

WGC.compute_rank(ignore_existing_files=True, sage="mod")
# WGC.plot_cohomology_dim(to_html=True)
# Euler char
WGC.print_dim_and_eulerchar()
WGC.print_cohomology_dim()


# WGC.plot_info()

# tt = WHairyGraphComplex.WHairyGraphVS(3,4,1,2)
# tt.build_basis(ignore_existing_files=True)
# tt.plot_all_graphs_to_file(skip_existing=False)
# tt.display_basis_plots()

# DSquareTestSingle(8,4,2,2)


# HG = WHairyGraphComplex.WHairyGraphVS(7,4,2,2)
# huntfor = "JSSqW@_?c??"
# huntforG = Graph(huntfor)

# ba = HG.get_basis_g6()
# print("In basis?: ", huntfor in ba)
# autom_list = huntforG.automorphism_group(partition=HG.get_partition()).gens()
# print("Odd automs?: ", HG._has_odd_automorphisms(huntforG,autom_list))
# print("auts", autom_list)
# print(huntforG)
# for p in autom_list:
#     print(list(p.tuple()))

# genlist=HG.get_generating_graphs()
# genlistg6 =[HG.graph_to_canon_g6(G)[0] for G in genlist]
# print("In gen list?: ", huntfor in genlistg6)
# testg6=HG.graph_to_canon_g6(huntforG)[0]
# print(testg6, testg6 in genlistg6)

# all_perm = [ list(range(0,HG.n_vertices+2)) + list(p) for p in itertools.permutations(range(HG.n_vertices+2, HG.n_vertices+HG.n_hairs+2)) ]
# print(all_perm)

# huntforG_h = huntforG.relabel({7:8, 8:7}, inplace=false)
# huntforG_h.add_edge(7,8)
# hairy_part= [list(range(HG.n_vertices+1)), list(range(HG.n_vertices+1, HG.n_vertices+HG.n_hairs+2))]
# huntforG_hg6 = huntforG_h.canonical_label(partition=hairy_part)

# somelist = HG._hairy_to_w_hair_2w(huntforG_h)
# somelistg6 =[HG.graph_to_canon_g6(G)[0] for G in somelist]
# print(huntfor in somelistg6)
# print(somelistg6)

# HG.build_basis(ignore_existing_files=True)

# check_graphs_vs_basis(tu.target, wwd)

# tt.domain.plot_all_graphs_to_file(skip_existing=False)
# tt.domain.display_basis_plots()
# tu.domain.plot_all_graphs_to_file(skip_existing=False)
# tu.domain.display_basis_plots()
# tu.target.plot_all_graphs_to_file(skip_existing=False)
# tu.target.display_basis_plots()


# g = Graph("EZq?")
# res = tt.operate_on(g)
# print(res)
# for G,v in res:
#     g6, s = tt.target.graph_to_canon_g6(G)
#     print(g6, v*s)


# class BasisTest(TestGraphComplex.BasisTest):
#     def setUp(self):
#         self.vs_list = [HairyGraphComplex.HairyGraphVS(v, l, h, even_edges, even_hairs) for (v, l, h, even_edges, even_hairs)
#                         in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


# class OperatorTest(TestGraphComplex.OperatorTest):
#     def setUp(self):
#         self.op_list = [HairyGraphComplex.ContractEdgesGO.generate_operator(v, l, h, even_edges, even_hairs) for
#                         (v, l, h, even_edges, even_hairs) in itertools.product(v_range, l_range, h_range, edges_types, hairs_types)]


# class GraphComplexTest(TestGraphComplex.GraphComplexTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]


# class CohomologyTest(TestGraphComplex.CohomologyTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]


# class SquareZeroTest(TestGraphComplex.SquareZeroTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, even_hairs, ['contract'])
#                         for (even_edges, even_hairs) in itertools.product(edges_types, hairs_types)]
#         self.gc_list += [HairyGraphComplex.HairyGC(v_range, l_range, h_range, even_edges, False, ['et1h']) for
#                          even_edges in edges_types]

# class AntiCommutativityTest(TestGraphComplex.AntiCommutativityTest):
#     def setUp(self):
#         self.gc_list = [HairyGraphComplex.HairyGC(v_range, l_range, h_range, False, False, ['contract', 'et1h'])]

# def suite():
#     suite = unittest.TestSuite()
#     suite.addTest(BasisTest('test_basis_functionality'))
#     suite.addTest(BasisTest('test_compare_ref_basis'))
#     suite.addTest(OperatorTest('test_operator_functionality'))
#     suite.addTest(OperatorTest('test_compare_ref_op_matrix'))
#     suite.addTest(GraphComplexTest('test_graph_complex_functionality'))
#     suite.addTest(CohomologyTest('test_cohomology_functionality'))
#     suite.addTest(SquareZeroTest('test_square_zero'))
#     suite.addTest(AntiCommutativityTest('test_anti_commutativity'))
#     return suite


# if __name__ == '__main__':
#     print("\n#######################################\n" + "----- Start test suite for hairy graph complex -----")
#     runner = unittest.TextTestRunner()
#     runner.run(suite())
