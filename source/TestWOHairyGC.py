# testing of basis-generation functionality using simple elementary examples and euler-characteristic
# for testing to be fully rigurous, the "data"-folder should be deleted before running the tests


import unittest
from sage.all import StandardTableaux
from WOHairyGC import WOHairyComponentGVS, WOHairyAggregatedGVS, WOHairyGVS, WOHairyGC


class TestBasisGeneration(unittest.TestCase):


    def test_double_leg_variants(self):
        # testing double-leg configurations
        
        WOHairyComponentGVS(n_vertices=0, n_loops=0, n=0, n_omega=2, n_epsilon=0).test_basis_len("omega-omega", basis_len=0)
        WOHairyComponentGVS(n_vertices=0, n_loops=0, n=0, n_omega=0, n_epsilon=2).test_basis_len("epsilon-epsilon", basis_len=1)
        WOHairyComponentGVS(n_vertices=0, n_loops=0, n=0, n_omega=1, n_epsilon=1).test_basis_len("omega-epsilon", basis_len=1)
        WOHairyComponentGVS(n_vertices=0, n_loops=0, n=1, n_omega=1, n_epsilon=0).test_basis_len("leg-omega", basis_len=1)
        WOHairyComponentGVS(n_vertices=0, n_loops=0, n=1, n_omega=0, n_epsilon=1).test_basis_len("leg-epsilon", basis_len=1)


        n_components = 1
        n_vertices = 0
        n_double_legs = 1

        WOHairyAggregatedGVS(genus=2, n=0, n_omega=2, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("omega-omega", basis_len=0)
        WOHairyAggregatedGVS(genus=2, n=0, n_omega=0, n_epsilon=2, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("epsilon-epsilon", basis_len=1)
        WOHairyAggregatedGVS(genus=2, n=0, n_omega=1, n_epsilon=1, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("omega-epsilon", basis_len=1)
        WOHairyAggregatedGVS(genus=1, n=1, n_omega=1, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("leg-omega", basis_len=1)
        WOHairyAggregatedGVS(genus=1, n=1, n_omega=0, n_epsilon=1, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("leg-epsilon", basis_len=1)


    def test_single_vertex_variants(self):
        # testing graphs with one vertex
        WOHairyComponentGVS(n_vertices=1, n_loops=0, n=0, n_omega=1, n_epsilon=2).test_basis_len("eps-omega-eps", basis_len=0)
        WOHairyComponentGVS(n_vertices=1, n_loops=0, n=0, n_omega=3, n_epsilon=0).test_basis_len("3-omega", basis_len=1)
        WOHairyComponentGVS(n_vertices=1, n_loops=0, n=0, n_omega=0, n_epsilon=3).test_basis_len("3-epsilon", basis_len=0)


        n_components = 1
        n_vertices = 1
        n_double_legs = 0

        WOHairyAggregatedGVS(genus=3, n=0, n_omega=1, n_epsilon=2, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("eps-omega-eps", basis_len=0)
        WOHairyAggregatedGVS(genus=3, n=0, n_omega=3, n_epsilon=0, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("3-omega", basis_len=1)
        WOHairyAggregatedGVS(genus=3, n=0, n_omega=0, n_epsilon=3, n_components=n_components, n_vertices=n_vertices, n_double_legs=n_double_legs).test_basis_len("3-epsilon", basis_len=0)


    def test_trees(self):
        # testing trees built from omegas and numbered legs on a single vertex

        n_components = 1
        n_vertices = 1
        n_epsilon = 0
        n_double_legs = 0

        for excess in range(6):

            for n_omega in range(excess + 4):
                
                genus = n_omega

                # excess = n_omega - 3 + 2*n -> 2n = excess - n_omega + 3
                two_n = excess - n_omega + 3

                if two_n % 2 == 0: 
                    n = int(two_n / 2)

                    # computing loop-order
                    n_hairs = n + n_omega + n_epsilon
                    n_edges = genus + n_vertices + n + n_double_legs - n_hairs - 1
                    n_loops = n_edges - n_vertices + (n_components - n_double_legs)

                    if n_omega + n >= 3 and n_omega >= 1:

                        WOHairyComponentGVS(n_vertices=n_vertices, 
                                            n_loops=n_loops, 
                                            n=n, n_omega=n_omega, n_epsilon=n_epsilon).test_basis_len("test_trees_WOHairyComponentGVS", basis_len=1)
                        
                        WOHairyAggregatedGVS(n_components=n_components, 
                                                    n_vertices=n_vertices, 
                                                    genus=genus, 
                                                    n=n, n_omega=n_omega, n_epsilon=n_epsilon,
                                                    n_double_legs=n_double_legs).test_basis_len("test_trees_WOHairyAggregatedGVS", basis_len=1)


    def test_mutliple_components(self):
        # testing graphs with multiple components ---

        WOHairyAggregatedGVS(genus=1, n=2, n_omega=2, n_epsilon=0, n_components=2, n_vertices=0, n_double_legs=2).test_basis_len("omega-1 & omega-2", basis_len=1)
        WOHairyAggregatedGVS(genus=5, n=2, n_omega=8, n_epsilon=0, n_components=4, n_vertices=2, n_double_legs=2).test_basis_len("omega-1 & omega-2 & 2x3-omega", basis_len=1)


    def test_bases_by_excess(self):

        # testing parameters (g, n, degree) in excess 0 with notrivial bases
        WOHairyGVS(genus=1, n=11, n_omega=11, degree=11).test_basis_len("excess 0: B_1_11", basis_len=1)
        WOHairyGVS(genus=3, n=8, n_omega=11, degree=14).test_basis_len("excess 0: B_3_8", basis_len=1)
        WOHairyGVS(genus=5, n=5, n_omega=11, degree=17).test_basis_len("excess 0: B_5_5", basis_len=1)
        WOHairyGVS(genus=7, n=2, n_omega=11, degree=20).test_basis_len("excess 0: B_7_2", basis_len=1)

        # excess 1
        WOHairyGVS(genus=2, n=10, n_omega=11, degree=13).test_basis_len("excess 1: B_2_10", basis_len=10)
        WOHairyGVS(genus=6, n=4, n_omega=11, degree=19).test_basis_len("excess 1: B_6_4", basis_len=5)

        # excess 2
        WOHairyGVS(genus=7, n=3, n_omega=11, degree=20).test_basis_len("excess 2: B_7_3", basis_len=16)



    @staticmethod
    def test_euler_char(genus, n, coefficients, diagrams, excess):

        print("testing for (g,n)=", (genus, n), "-----------------------")

        computed_excess = 3*genus + 2*n - 25

        if computed_excess < 0: assert excess < 0
        else: assert computed_excess == excess

        euler_char = WOHairyGVS.compute_euler_char(genus, n)
        
        assert len(coefficients) == len(diagrams)

        euler_char_reference = 0
        for coefficient, diagram in zip(coefficients, diagrams):
            euler_char_reference += coefficient * StandardTableaux(diagram).cardinality()

        print("euler_char:", euler_char)
        print("euler_char_reference:", euler_char_reference)
        assert euler_char == euler_char_reference, (euler_char, euler_char_reference)


    def test_eulerChar(self):
        # excess < 0
        for i in range(5): self.test_euler_char(5, i, coefficients=[], diagrams=[], excess=-1)
        for i in range(4): self.test_euler_char(6, i, coefficients=[], diagrams=[], excess=-1)
        for i in range(2): self.test_euler_char(7, i, coefficients=[], diagrams=[], excess=-1)
        self.test_euler_char(8, 0, coefficients=[], diagrams=[], excess=-1)

        # excess = 0
        self.test_euler_char(5, 5, coefficients=[-1], diagrams=[[1,1,1,1,1]], excess=0)
        self.test_euler_char(7, 2, coefficients=[1], diagrams=[[1,1]], excess=0)

        # excess = 1
        self.test_euler_char(6, 4, coefficients=[-1], diagrams=[[2,1,1]], excess=1)
        self.test_euler_char(8, 1, coefficients=[], diagrams=[], excess=1)

        # excess = 2
        self.test_euler_char(5, 6, coefficients=[-1, 2], diagrams=[[1,1,1,1,1,1], [3,1,1,1]], excess=2)
        self.test_euler_char(7, 3, coefficients=[1, -2], diagrams=[[1,1,1], [3]], excess=2)
        self.test_euler_char(9, 0, coefficients=[1], diagrams=[[]], excess=2)

        # excess = 3
        self.test_euler_char(8, 2, coefficients=[1], diagrams=[[2]], excess=3)
        self.test_euler_char(6, 5, coefficients=[-1, 1, 1, 2], diagrams=[[2,1,1,1], [3,1,1], [3,2], [4,1]], excess=3)

        # excess = 4
        self.test_euler_char(9, 1, coefficients=[1], diagrams=[[1]], excess=4)
        self.test_euler_char(7, 4, coefficients=[2, -2], diagrams=[[2,2], [3,1]], excess=4)


        # excess = 5
        self.test_euler_char(10, 0, coefficients=[-2], diagrams=[[]], excess=5)
        self.test_euler_char(8, 3, coefficients=[4], diagrams=[[1,1,1]], excess=5)
        #self.test_euler_char(6, 6, coefficients=[-4, -2, -3, -3, 1, -2, -2, 2, -3, -2 ,-3], 
        #                diagrams=[[1,1,1,1,1,1], [2,1,1,1,1], [2,2,1,1], [2,2,2], 
        #                        [3,1,1,1], [3,2,1], [3,3], [4,1,1], [4,2], [5,1], [6]], excess=5)

        # excess = 6
        self.test_euler_char(9, 2, coefficients=[4], diagrams=[[2]], excess=6)
        #self.test_euler_char(7, 5, coefficients=[-1, -7, -6, 2, 1, 4], diagrams=[[1,1,1,1,1], [2,1,1,1], [3,1,1], [3,2], [4,1], [5]], excess=6)

        # excess = 7
        self.test_euler_char(10, 1, coefficients=[-6], diagrams=[[1]], excess=7)
        #self.test_euler_char(8, 4, coefficients=[10, 2, -2, -10, -3], diagrams=[[1,1,1,1], [2,1,1], [2,2], [3,1], [4]], excess=7)

        # excess = 8
        self.test_euler_char(11, 0, coefficients=[2], diagrams=[[]], excess=8)
        #self.test_euler_char(9, 3, coefficients=[1, 11, 4], diagrams=[[1,1,1], [2,1], [3]], excess=8)

        # excess = 9
        #self.test_euler_char(10, 2, coefficients=[-9, -2], diagrams=[[1,1], [2]], excess=9)

        # excess = 10
        #self.test_euler_char(11, 1, coefficients=[1], diagrams=[[1]], excess=10)




class TestOperators(unittest.TestCase):

    

    def DSquareTest_EpsToOmega(self):
        g_n_pairs = [
            (7, 2), (5, 5), (3, 8), (1, 11),
            (8, 1), (6, 4), (4, 7), (2, 10),
            (9, 0), (7, 3), (5, 6), (3, 9), (1, 12),
            (8, 2), (6, 5), (4, 8), (2, 11),
            (9, 1), (7, 4), (5, 7), (3, 10),(1, 13),
            (10, 0), (8, 3), (9, 2), (10, 1), (11, 0)
            ]
        for (genus, n) in g_n_pairs:
            WOHairyGC.DSquareTest_single(operator='epstoomega', genus=genus, n=n)


    def DSquareTest_ContractEdges(self):
        g_n_pairs = [
            (7, 2), (5, 5), (3, 8), (1, 11),
            (8, 1), (6, 4), (4, 7), (2, 10),
            (9, 0), (7, 3), (5, 6), (3, 9), (1, 12),
            (8, 2), (6, 5), (4, 8), (2, 11)
            #,(9, 1), (7, 4), (5, 7), (3, 10),(1, 13)
            #,(10, 0), (8, 3), (9, 2), (10, 1), (11, 0)
            ]
                
        for (genus, n) in g_n_pairs:
            WOHairyGC.DSquareTest_single(operator='contract', genus=genus, n=n)


    def Anticommutativity_Test(self):
        g_n_pairs = [
            (7, 2), (5, 5), (3, 8), (1, 11),
            (8, 1), (6, 4), (4, 7), (2, 10),
            (9, 0), (7, 3), (5, 6), (3, 9), (1, 12),
            (8, 2), (6, 5), (4, 8), (2, 11)
            #,(9, 1), (7, 4), (5, 7), (3, 10),(1, 13)
            #,(10, 0), (8, 3), (9, 2), (10, 1), (11, 0)
            ]
        for (genus, n) in g_n_pairs:
            WOHairyGC.Anticomm_Test_single(genus=genus, n=n)


    @staticmethod
    def cohomology_dim_test_single(genus, n, nonzero_degrees, diagram_lists, n_omega = 11):

        assert len(nonzero_degrees) == len(diagram_lists)

        print("(genus, n) = ", (genus, n))

        deg_min = 22 - n_omega + genus - 1
        deg_max = 3*genus + n + 19 - 2*n_omega
        degree_range = range(deg_min, deg_max+3)
        
        r2 = None
        for degree in degree_range:
            
            cohom_dim, r2 = WOHairyGC.compute_cohomology_dim(degree=degree, genus=genus, n=n, n_omega=n_omega, prev_r2=r2)
            
            if cohom_dim == 0:
                assert degree not in nonzero_degrees, "(g,n) = "+ str((genus, n)) + " degree: " + str(degree)
                pass
            else:
                assert cohom_dim > 0
                assert degree in nonzero_degrees, "(g,n) = "+ str((genus, n)) + " degree: " + str(degree)
                index = nonzero_degrees.index(degree)
                diagram_list = diagram_lists[index]

                test_dim = 0
                for digaram in diagram_list:
                    test_dim += StandardTableaux(digaram).cardinality()
                
                assert cohom_dim == test_dim, "(g,n) = "+ str((genus, n))



    def Cohom_dim_Test(self):
        g_n_dim_pairs = [   

            # excess 0
            (7, 2, [20], [[[1,1]]]),
            (5, 5, [17], [[[1,1,1,1,1]]]),
            (3, 8, [14], [[[1,1,1,1,1,1,1,1]]]),
            (1, 11, [11], [[[1,1,1,1,1,1,1,1,1,1,1]]]),


            # excess 1
            (8, 1, [], []),
            (6, 4, [19], [[[2,1,1]]]),
            (4, 7, [16], [[[2,1,1,1,1,1]]]),
            (2, 10, [13], [[[2,1,1,1,1,1,1,1,1]]]),


            # excess 2
            (9, 0, [22], [[[]]]),
            (7, 3, [20, 21], [[[1,1,1]], [[3], [3]]]),
            (5, 6, [17, 18], [[[1,1,1,1,1,1]], [[3,1,1,1], [3,1,1,1]]]),
            (3, 9, [14, 15], [[[1,1,1,1,1,1,1,1,1]], [[3,1,1,1,1,1,1], [3,1,1,1,1,1,1]]]),
            (1, 12, [12], [[[3,1,1,1,1,1,1,1,1,1]]]),
            

            # excess 3
            (8, 2, [22], [[[2]]]),
            (6, 5, [19, 20], [[[2,1,1,1]], 
                              [[4,1], [4,1], [3,2], [3,1,1]]]),
            (4, 8, [16, 17], [[[2,1,1,1,1,1,1]], 
                              [[4,1,1,1,1], [4,1,1,1,1], [3,2,1,1,1], [3,1,1,1,1,1]]]),
            (2, 11, [14], [[[4,1,1,1,1,1,1,1], [3,2,1,1,1,1,1,1], [3,1,1,1,1,1,1,1,1]]]),


            # excess 4
            (9, 1, [22], [[[1]]]), # TODO: Mistake in Paper: "k=24 instead of k=22" ?
            (7, 4, [21, 22], [[[3,1], [3,1]], 
                              [[2,2], [2,2]]]),
            (5, 7, [18, 19], [[[3,1,1,1,1], [3,1,1,1,1]], 
                              [[5,1,1], [5,1,1], [5,1,1], [4,2,1], [4,1,1,1], [3,3,1], [3,3,1], [3,2,1,1], [3,2,1,1], [2,2,1,1,1], [2,2,1,1,1]]]),
            (3, 10, [15, 16], [[[3,1,1,1,1,1,1,1]], 
                              [[5,1,1,1,1,1], [5,1,1,1,1,1], [4,2,1,1,1,1], [4,1,1,1,1,1,1], [3,3,1,1,1,1], [3,3,1,1,1,1], [3,2,1,1,1,1,1], [3,2,1,1,1,1,1], [2,2,1,1,1,1,1,1], [2,2,1,1,1,1,1,1]]]),
            (1, 13, [13], [[[5,1,1,1,1,1,1,1,1], [3,3,1,1,1,1,1,1,1], [3,2,1,1,1,1,1,1,1,1], [2,2,1,1,1,1,1,1,1,1,1]]]),

        ]
        
        for (genus, n, nonzero_degrees, diagram_lists) in g_n_dim_pairs:
            TestOperators.cohomology_dim_test_single(genus=genus, n=n, nonzero_degrees=nonzero_degrees, diagram_lists=diagram_lists)


        
    


def suite():

    suite = unittest.TestSuite()

    suite.addTest(TestBasisGeneration('test_double_leg_variants'))
    suite.addTest(TestBasisGeneration('test_single_vertex_variants'))
    suite.addTest(TestBasisGeneration('test_trees'))
    suite.addTest(TestBasisGeneration('test_mutliple_components'))
    suite.addTest(TestBasisGeneration('test_bases_by_excess'))

    suite.addTest(TestBasisGeneration('test_eulerChar'))


    suite.addTest(TestOperators('DSquareTest_EpsToOmega')) 
    suite.addTest(TestOperators('DSquareTest_ContractEdges')) 
    suite.addTest(TestOperators('Anticommutativity_Test'))

    suite.addTest(TestOperators('Cohom_dim_Test'))

    return suite



if __name__ == '__main__':
    print("\n#######################################\n" + "----- Start test suite WOHairy Baisis Generation -----")
    runner = unittest.TextTestRunner(verbosity=2)
    runner.run(suite())
