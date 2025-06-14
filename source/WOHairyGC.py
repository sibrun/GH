"""
Author: Pascal Skipness

We consider Graph complexes based on simple graphs with numbered hairs and hairs of two (omega- and epsilon-)decorations
as in the graph complex computing weight 11 cohomology. The omega decorations are considered odd.
For more detailed Information see the following references:
- Weight 11 compactly supported cohomology of moduli spaces of curves: https://arxiv.org/abs/2302.04204
- Weight 11 Compactly Supported Cohomology of Moduli Spaces of Curves in excess four: https://arxiv.org/abs/2407.16372

This file implements the basis-generation and differential operators of these graph complexes.
This is achieved using the following classes:

    WOHairyComponentGVS: 
        generates the connected compontents of such graphs.
        This is achieved as follows:
        1. generate Hairy graphs with numbered hairs using CHairyGraphComplex.py
        2. relable some of the hairs with omega or epsilon

    WOHairyAggregatedGVS:    
        recursively joins together connected graphs from WOHairyComponentGVS 
        to get basis-graphs with any number of connected components.

    WOHairyGVS:
        selects, for some given cohomological degree, the relevant graphs 
        from WOHairyAggregatedGVS and puts them in a single file.

    EpsToOmegaGO:
        graph operator which replaces an epsilon-hair with an omega-hair
        main functionality implemented in operate_on

    ContractEdgeGO:
        graph operator which contracts an edge of the graph
        main functionality implemented in operate_on

    WOHairyGC:
        graph complex which uses the above classes to generate the basis and the differential operators
        further contains testing functions and cohomology computation functions

Unittests are implemented in the file "TestWOHairyGC.py"

"""



from sage.all import *
import Shared
import Parameters
import CHairyGraphComplex
import itertools
import math
import GraphVectorSpace
import matplotlib.pyplot as plt
import SymmetricGraphComplex
import GraphOperator
import os
import Shared
from copy import copy
import GraphComplex
import numpy as np
import csv

# ------- Graph Vector Space --------
class WOHairyComponentGVS(CHairyGraphComplex.CHairyGraphVS):

    def __init__(self, n_vertices, n_loops, n, n_omega, n_epsilon):

        self.n = n
        self.n_vertices = n_vertices
        self.n_loops = n_loops # inner loops! (without merging epsilons and omegas)
        self.n_omega = n_omega
        self.n_epsilon = n_epsilon
        self.even_edges = False
        self.n_hairs = self.n + self.n_omega + self.n_epsilon

        # computation of inner-edges
        if self.n_vertices > 0: self.n_edges = self.n_loops + self.n_vertices - 1 # case of no double-leg
        else: self.n_edges = 0 # case of double-leg


    def __hash__(self):
        return hash("wo_comp_gra%d_%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wo_comp_gra%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, "wohairy_component", s)

    def get_ref_basis_file_path(self):
        s = "wo_comp_gra%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, "wohairy_component", self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('vertices', self.n_vertices), ('loops', self.n_loops), ('numbered_hairs', self.n), 
                                   ('omegas', self.n_omega), ('epsilons', self.n_epsilon)])


    def get_partition(self):

        if self.n_vertices > 0: inner_vertices = [list(range(0, self.n_vertices))]
        else: inner_vertices = []

        numbered_vertices = [[j] for j in range(self.n_vertices, self.n_vertices + self.n)]

        omega_vertices = [list(range(self.n_vertices + self.n, self.n_vertices + self.n + self.n_omega))]

        if self.n_epsilon > 0: epsilon_vertices = [list(range(self.n_vertices + self.n + self.n_omega, 
                                                            self.n_vertices + self.n + self.n_omega + self.n_epsilon))]
        else: epsilon_vertices = []
        
        partition = inner_vertices + numbered_vertices + omega_vertices + epsilon_vertices

        return partition



    def is_valid(self):

        # allowed values
        l = self.n_vertices >= 0 and self.n_loops >= 0 and self.n_edges >= 0 \
        and self.n >= 0 and self.n_omega >= 0 and self.n_epsilon >= 0
        # At least trivalent internal vertices.
        l = l and (3 * self.n_vertices <= 2 * self.n_edges + self.n_hairs)
        # At most a full graph.
        l = l and self.n_edges <= (self.n_vertices) * (self.n_vertices - 1) / 2
        # connected graph must contain at least one omega or epsilon
        l = l and self.n_omega >= 1 or self.n_epsilon >= 1
        # if double-leg
        if self.n_vertices == 0: l = l and self.n_hairs == 2

        return l
    

    def get_generating_graphs(self):

        if not self.is_valid():
            print("self is not valid")
            return []
    
        # produce all hairy graphs
        hairy_graphs = self.get_hairy_graphs(self.n_vertices, self.n_loops, self.n_hairs, include_novertgraph=True)

        # produce all neccesary permutations of the hairs
        # for more information look at the function "multiset_permutations" below
        all_perm = multiset_permutations(n_vertices=self.n_vertices, n=self.n, n_omega=self.n_omega, n_epsilon=self.n_epsilon)

        return (G.relabel(p, inplace=False) for G in hairy_graphs for p in all_perm)


    
    def perm_sign(self, G, p):
        
        #inputs:
        #    G: sage-graph with vertices labled [0, ..., n-1] 
        #    p: permutation of the vertices given by a list [p(0), ..., p(n-1)]
        #output: sign of the permutation induced on the edges of the Graph G
        
        assert set(G.vertices()) == set(p)
        
        #print("vertex-permutation: ", p)

        G1 = copy(G)
        Shared.enumerate_edges(G1)
        # We permute the graph, and read of the new labels
        G1.relabel(p, inplace=True)
        sgn = Shared.Perm([j for (u, v, j) in G1.edges(sort=True)]).signature()
        
        # Compute the extra contribution from omega-hairs.
        if self.n_omega > 0:
            omega_hairs = p[self.n_vertices + self.n : self.n_vertices + self.n + self.n_omega]
            #print("omega_hairs: ", omega_hairs)
            sgn_omega_perm = Shared.Perm.shifted(omega_hairs).signature()
            #print("sign of induced omega permutation: ", sgn_omega_perm)
            sgn *= sgn_omega_perm

        return sgn

    def test_basis_len(self, test_name, basis_len):

        if basis_len > 0: assert self.is_valid(), test_name

        self.build_basis(ignore_existing_files=True)

        assert self.get_dimension() == basis_len, test_name
    

# helper-functions ---

def multiset_permutations(n_vertices, n, n_omega, n_epsilon):
    """
    this is a helper-function for WOHairyComponentGVS.get_generating_graphs()
    instead of considering all permutations of the hairs of the graph, we can make the following restriction:
    - omega-hairs are identical and hence need not be permuted with each other
    - the same holds for the epsilon-hairs

    Using this we need only to consider (n + n_omega + n_epsilon)! / (n_omega! * n_epsilon!)
    instead of (n + n_omega + n_epsilon)! permutations
    
    :Example:

    - input: [n_vertices=0, n=1, n_omega=2, n_epsilon=0]
    - possible permutations of the hairs: [1, omega, omega], [omega, 1, omega], [omega, omega, 1]
    - intuitively the output should be: [0, 1, 2], [1, 0, 2], [2, 0, 1]
    - BUT: the partition is always given by [[0], [1,2]]
      hence the new "1"-vertex needs to get relabled to 0
    - because of this the output consists of the inverted permutations of the intuitive output:
      [0, 1, 2], [1, 0, 2], [1, 2, 0]
    """
    permutations = []

    L = n + n_omega + n_epsilon
    
    # 1) Choose positions for the n_omega identical omega's
    for omega_positions in itertools.combinations(range(L), n_omega):
        set_omega = set(omega_positions)
        
        # 2) Choose positions (from the leftover) for the n_epsilon identical epsilon's
        leftover_after_omega = [i for i in range(L) if i not in set_omega]
        for epsilon_positions in itertools.combinations(leftover_after_omega, n_epsilon):
            set_epsilon = set(epsilon_positions)
            
            # 3) The remaining L - n_omega - n_epsilon positions are for the n distinct items.
            leftover_for_distinct = [i for i in range(L)
                                     if i not in set_omega and i not in set_epsilon]
            
            # Permute the n distinct items in all possible ways
            for perm in itertools.permutations(leftover_for_distinct):
                
                # compute the full permutation
                full_perm = list(perm) + list(omega_positions) + list(epsilon_positions) 
                assert set(full_perm) == set(list(range(L)))

                # invert the permutation
                full_perm_inv = [None]*L
                for i in range(L):
                    full_perm_inv[full_perm[i]] = i
                
                # add the identity-permutaion for the inner vertices and
                # translate the permutation by the number of inner vertices 
                full_perm_inv_translated = list(range(0, n_vertices)) + [n_vertices + i for i in full_perm_inv] 

                permutations.append(full_perm_inv_translated)

    # assert that we have generated the correct number of permutations
    assert len(permutations) == math.factorial(n + n_omega + n_epsilon) / (math.factorial(n_omega) * math.factorial(n_epsilon))

    return permutations











class WOHairyAggregatedGVS(WOHairyComponentGVS):

    def __init__(self, n_components, n_vertices, genus, n, n_omega, n_epsilon, n_double_legs):

        self.n_components = n_components
        self.n = n
        self.n_vertices = n_vertices
        self.genus = genus
        self.n_omega = n_omega
        self.n_epsilon = n_epsilon
        self.n_double_legs = n_double_legs
        self.n_hairs = self.n + self.n_omega + self.n_epsilon
        self.total_num_vertices = self.n_vertices + self.n_hairs
        self.excess = 3*(self.genus - 1) + 2*self.n - 2*self.n_omega

        # computation of the inner edges by genus
        # g = (E-V+num_components) + 1 of the contracted graph:
        # g = (n_edges + n_hairs - n_double_legs) - (n_vertices + n + 1) + 1 + 1 
        #   = n_edges + n_hairs - n_vertices - n - n_double_legs + num_components
        self.n_edges = genus + n_vertices + n + n_double_legs - self.n_hairs - 1

        # loop-order of the inner graph
        self.n_loops = self.n_edges - n_vertices + (self.n_components - self.n_double_legs)



    def cohom_degree(self):
        # deg = 22 + (all_edges - n_omega - n)
        # where all_edges = inner_edges + n + n_eps + n_omega - n_double_legs
        # -> deg = 22 + inner_edges + n_eps - n_double_legs
        return 22 + self.n_edges + self.n_epsilon - self.n_double_legs
    

    def __hash__(self):
        return hash("wo_agg_gra%d_%d_%d_%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wo_agg_gra%d_%d_%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, "wohairy_aggregated", s)

    def get_ref_basis_file_path(self):
        s = "wo_agg_gra%d_%d_%d_%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, "wohairy_aggregated", self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('components', self.n_components), ('vertices', self.n_vertices), ('genus', self.genus), ('numbered_hairs', self.n), 
                                   ('omegas', self.n_omega), ('epsilons', self.n_epsilon), ('double_legs', self.n_double_legs)])


    def is_valid(self):
        # possible values
        l = self.n_components >= 1 \
        and self.n_vertices >= 0 \
        and self.genus >= 1 \
        and self.n >= 0 \
        and self.n_omega >= 0 \
        and self.n_epsilon >= 0 \
        and self.n_double_legs >= 0 \
        and self.n_hairs >= 0 \
        and self.n_edges >= 0 \
        and self.n_loops >= 0 \
        and self.excess >= 0
        #print("is_valid() after value-checking:", l)
        # each n_double_leg corresponds to a connected component 
        l = l and self.n_double_legs <= self.n_components
        # number of components which are not double-legs: self.n_components - self.n_double_legs
        l = l and self.n_vertices >= self.n_components - self.n_double_legs
        # If all components are double-legs:
        if self.n_components == self.n_double_legs: 
            l = l and self.n_vertices == 0
            l = l and self.n_hairs == 2*self.n_double_legs
            l = l and self.n_loops == 0
            l = l and self.n_edges == 0
        
        # each connected component must at least contain one omega or epsilon
        l = l and self.n_omega + self.n_epsilon >= self.n_components
        # At least trivalent internal vertices. (each double_leg removes two hair-decorations from the equation)
        l = l and (3 * self.n_vertices <= 2 * self.n_edges + (self.n_hairs - 2*self.n_double_legs)) 
        # At most a full graph.
        l = l and self.n_edges <= (self.n_vertices) * (self.n_vertices - 1) / 2
        # each double-leg contains exactly two hair lables (which are not inner vertices)
        l = l and self.n_double_legs <= (self.total_num_vertices - self.n_vertices) / 2
        # each other connected component contains at least one vertex
        l = l and self.n_components - self.n_double_legs <= self.total_num_vertices - 2*self.n_double_legs
        # due to symmetry reasons, each vertex can have maximally one epsilon-hair attached to it otherwise the graph is equal to zero
        l = l and self.n_epsilon <= 2*self.n_double_legs + self.n_vertices
        # cannot have more eps than excess
        l = l and self.n_epsilon <= self.excess
        
        # test for single-vertex trees with omega and numbered hairs
        if self.n_vertices == 1 and self.n_components == 1 and self.n_epsilon == 0:
            l = l and self.n_double_legs == 0
            l = l and (self.excess - self.n_omega + 3) % 2 == 0
            l = l and self.excess == self.n_omega - 3 + 2*self.n
            l = l and self.n_omega + self.n >= 3
        

        # lemma: for excess <= 4, the loop order is 0
        if self.excess <= 4: l = l and self.n_loops == 0
        
        return l


    

    def get_generating_graphs(self):
        #print("Aggregated: get_generating_graphs")
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            print("self is not valid")
            return []
    

        # recursive aggregation by n_connected components
        assert isinstance(self.n_components, int)
        assert (self.n_components >= 1)

        # if we only have one connected component the definition is equal to the definition of ModifedCHairyGraphComplex
        if self.n_components == 1:

            V = WOHairyComponentGVS(n_vertices=self.n_vertices,
                                                              n_loops=self.n_loops,
                                                              n=self.n,
                                                              n_omega=self.n_omega,
                                                              n_epsilon=self.n_epsilon)

            V.build_basis(ignore_existing_files=False)

            #print("prebuilt basis len:", len([G for G in V.get_basis()]))

            for G in V.get_basis():

                assert len(G.vertices()) == self.total_num_vertices
                assert len(G.edges()) == self.n_edges + self.n_hairs - self.n_double_legs

                yield G


        # else we need to iterate through different possible configurations of joining 
        # a graph complex with n_components-1 connected components with one with only a single connected component.
        else: 
            #print("recursion-step for n_components=" + str(self.n_components))

            configurations = self.get_configurations()
            #print(configurations)
            #print("number of combinatorically possible configurations: " + str(len(configurations)))
            #print(configurations)
            for configuration in configurations:
                #print("configuration: " + str(configuration))
                (n_vert_1, n_vert_2, genus_1, genus_2, n_1, n_2, n_omega_1, n_omega_2, n_eps_1, n_eps_2, n_DL_1, n_DL_2) = configuration

                n_components_1 = 1
                n_components_2 = self.n_components - 1
                
                V_1 = WOHairyAggregatedGVS(n_components=n_components_1,
                                         n_vertices=n_vert_1,
                                         genus=genus_1,
                                         n=n_1,
                                         n_omega=n_omega_1,
                                         n_epsilon=n_eps_1,
                                         n_double_legs=n_DL_1)
                
                if not V_1.is_valid(): continue

                V_2 = WOHairyAggregatedGVS(n_components=n_components_2,
                                         n_vertices=n_vert_2,
                                         genus=genus_2,
                                         n=n_2,
                                         n_omega=n_omega_2,
                                         n_epsilon=n_eps_2,
                                         n_double_legs=n_DL_2)
                
                if not V_2.is_valid(): continue

                # check if excess adds up
                if self.excess != V_1.excess + V_2.excess: continue

                V_1.build_basis(ignore_existing_files=False)
                #print("built basis 1")
                V_2.build_basis(ignore_existing_files=False)
                #print("built basis 2")
                #print("dimension V_1:", V_1.get_dimension())
                #print("dimension V_2:", V_2.get_dimension())
                
                # join graphs
                #print("joining possible graph combinations")
                for G_single_comp in V_1.get_basis():
                    for G2_mult_comp in V_2.get_basis():
                        #print("joining graphs")

                        # join graphs:
                        # G now has vertices (0,i) and (1,j)
                        G = G_single_comp.disjoint_union(G2_mult_comp)

                        assert len(G.vertices()) == self.total_num_vertices
                        assert len(G.edges()) == self.n_edges + self.n_hairs - self.n_double_legs

                        # relabel vertices to 0,1,2,...
                        G.relabel({old: new for new, old in enumerate(G.vertices())}, inplace=True)
                
                        # permute vertices accordingly to how the partition is made
                        # i.e. we need to rearrange the vertices according to the partition
                        # AND we need to permute the numbered hairs between the two components
                        # this is done by the functions reorder_vertices and get_cross_permutations
                        old_order = G.vertices()
                        new_orders = reorder_vertices(old_order, n_vert_1, n_vert_2, n_1, n_2, n_omega_1, n_omega_2, n_eps_1, n_eps_2)

                        for new_order in new_orders:
                            mapping = {new:old for old, new in zip(old_order, new_order)}
                            
                            G_new = G.copy()
                            G_new.relabel(mapping)

                            yield G_new

                            """
                            print("G.vertices(): " + str(G.vertices()))
                            print("G.edges(): " + str(G.edges()))
                            print("H.vertices(): " + str(H.vertices()))
                            print("H.edges(): " + str(H.edges()))
                            """
        
    


    def get_configurations(self):
        
        configurations = []
        
        for n_vert_1 in range(0, self.n_vertices + 1):
            n_vert_2 = self.n_vertices-n_vert_1
            #print("n_vert_1: " + str(n_vert_1))

            # ensure that genus_1 and genus_2 are both >= 1
            for genus_1 in range(1, self.genus+1):
                # g = (E-V+1) + 1 of the contracted graph:
                #   = (n_edges + n_hairs - n_double_legs) - (n_vertices + n + 1) + 1 + 1 
                #   = n_edges + n_hairs - n_vertices - n_double_legs - n + 1 
                # g_tot = g_1 + g_2 - 1
                genus_2 = self.genus - genus_1 + 1
                assert genus_1 >= 1 and genus_2 >= 1
                #print("genus_1: " + str(genus_1))

                for n_1 in range(self.n + 1):
                    n_2 = self.n-n_1
                    #print("n_1: " + str(n_1))

                    for n_omega_1 in range(self.n_omega + 1):
                        n_omega_2 = self.n_omega-n_omega_1
                        #print("n_omega_1: " + str(n_omega_1))

                        for n_epsilon_1 in range(self.n_epsilon + 1):
                            n_epsilon_2 = self.n_epsilon-n_epsilon_1
                            #print("n_epsilon_1: " + str(n_epsilon_1))

                            for n_double_legs_1 in range(0, self.n_double_legs+1):
                                n_double_legs_2 = self.n_double_legs - n_double_legs_1
                            
                                if (n_omega_1 >= 1 or n_epsilon_1 >= 1) and (n_omega_2 >= 1 or n_epsilon_2 >= 1):

                                    configurations.append((n_vert_1, n_vert_2,
                                                        genus_1, genus_2,
                                                        n_1, n_2,
                                                        n_omega_1, n_omega_2,
                                                        n_epsilon_1, n_epsilon_2,
                                                        n_double_legs_1, n_double_legs_2))
        
        return configurations 

        

# helper-functions ---
        
def list_elems_from_range(input_list, num, dist):
    return [input_list[i] for i in list(range(num, num + dist))]

def list_elems_from_double_range(input_list, num1, dist1, num2, dist2):
    """
    input: some list 

    output:

    two disjoint slices of the list (starting at num_i and of length dist_i)
    merged together to a single new list
    """
    assert max(num1 + dist1, num2 + dist2) <= len(input_list)
    assert num1 + dist1 <= num2 or num2 + dist2 <= num1 # disjointness of the slice

    return list_elems_from_range(input_list, num1, dist1) + list_elems_from_range(input_list, num2, dist2)


def get_cross_permutations(list_1, list_2):
    """
    Since the individual components already contain all of the internal permutations of the numbered hairs, 
    we only need to consider order-preserving permutations which swap numbered hairs between the two components.
    Let n_1, n_2 be the respective number of numbered hairs
    Then within the generation of the individual components we have already considered 
    the n_i! permutations of the elements in list_i for i=1,2.
    Hence in order to consider all (n_1+n_2)! permutations it is only left to check 
    additional (n_1+n_2)! / (n_1!*n_2!) permutations for any combination of such components.

    More operationally we compute the following:
    Generate all permutations of S that satisfy the constraints:

    - S is partitioned into disjoint subsets A and B.
    - Each element in A is either fixed or swapped with an element in B (and vice versa).
    - If two elements i, j in A (i < j) are swapped, their images in B preserve order (σ(i) < σ(j)).

    Likewise, if two in B are swapped, their images in A preserve order.

    Yields each valid permutation as a tuple of length n (in one-line notation).
    """
        # assert that list_1 and list_2 are disjoint with n combined elements
    assert len(list_1) + len(list_2) == len(set(list_1 + list_2))

    permutations = []

    # Iterate over possible swap counts
    for k in range(0, min(len(list_1), len(list_2)) + 1):
        # Choose k elements from A and k from B to swap
        for combA in itertools.combinations(list_1, k):
            #print("combA", combA)
            for combB in itertools.combinations(list_2, k):
                #print("combB", combB)
                # Pair the chosen elements in sorted order (combinations are already sorted)
                perm = list_1 + list_2  # start with identity permutation
                for a_elem, b_elem in zip(combA, combB):
                    # Swap a_elem and b_elem in the permutation
                    #print("before swap", perm)
                    i = perm.index(a_elem)
                    #print(i)
                    j = perm.index(b_elem)
                    #print(j)
                    
                    # Swap in place
                    perm[i], perm[j] = perm[j], perm[i]
                    #print("after swap", perm)
                    
                permutations.append(perm)
    
    
    n_1 = len(list_1)
    n_2 = len(list_2)
    assert len(permutations) == math.factorial(n_1 + n_2) / (math.factorial(n_1) * math.factorial(n_2))

    return permutations


def reorder_vertices(old_order, n_vert_1, n_vert_2, n_1, n_2, n_omega_1, n_omega_2, n_epsilon_1, n_epsilon_2):
    assert len(old_order) == n_vert_1 + n_vert_2 + n_1 + n_2 + n_omega_1 + n_omega_2 + n_epsilon_1 + n_epsilon_2

    new_orders = []
    #print("old_order: " + str(old_order))
    divider = n_vert_1 + n_1 + n_omega_1 + n_epsilon_1

    inner_vertices = list_elems_from_double_range(old_order, 0, n_vert_1, divider, n_vert_2)
    #print("inner_vertices: " + str(inner_vertices))
    
    numbered_vertices_1 = list_elems_from_range(old_order, n_vert_1, n_1)
    numbered_vertices_2 = list_elems_from_range(old_order, divider+n_vert_2, n_2)
    #print("numbered_vertices: " + str(numbered_vertices))

    omega_vertices = list_elems_from_double_range(old_order, n_vert_1+n_1, n_omega_1, divider+n_vert_2+n_2, n_omega_2)
    #print("omega_vertices: " + str(omega_vertices))

    epsilon_vertices = list_elems_from_double_range(old_order, n_vert_1+n_1+n_omega_1, n_epsilon_1, divider+n_vert_2+n_2+n_omega_2, n_epsilon_2)
    #print("epsilon_vertices: " + str(epsilon_vertices))

    #print("list1", numbered_vertices_1)
    #print("list2", numbered_vertices_2)
    #for permuted_numbered_verticies in itertools.permutations(numbered_vertices_1 + numbered_vertices_2):
    #print("un-permuted:", numbered_vertices_1 + numbered_vertices_2)
    for permuted_numbered_verticies in get_cross_permutations(numbered_vertices_1, numbered_vertices_2):
        
        #print("permuted_numbered_verticies: " + str(permuted_numbered_verticies))
        new_order = inner_vertices + list(permuted_numbered_verticies) + omega_vertices + epsilon_vertices

        #print("new_order: " + str(new_order))
        assert len(new_order) == len(old_order)
        assert set(new_order) == set(range(len(new_order)))

        new_orders.append(new_order)

    return new_orders






class WOHairyGVS(WOHairyAggregatedGVS):

    def __init__(self, genus, n, n_omega, degree):

        self.genus = genus
        self.n = n
        self.n_omega = n_omega
        self.degree = degree   
        
        self.n_vertices = degree - 22 + n_omega - genus + 1
        self.excess = 3*genus + 2*n - 25 

    def get_type(self):
        return 'wohairy'
    
    def __eq__(self, other):
        return self.genus == other.genus \
                and self.n == other.n \
                and self.n_omega == other.n_omega \
                and self.degree == other.degree

    def get_n_epsilon_from_graph(self, G):

        n_epsilon = len(G.vertices()) - self.n_vertices - self.n - self.n_omega

        assert isinstance(n_epsilon, int)
        assert n_epsilon >= 0

        return n_epsilon


    def __hash__(self):
        return hash("wo_fin_gra%d_%d_%d_%d" % self.get_ordered_param_dict().get_value_tuple())

    def get_basis_file_path(self):
        s = "wo_fin_gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, "wohairy", s)

    def get_ref_basis_file_path(self):
        s = "wo_fin_gra%d_%d_%d_%d.g6" % self.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, "wohairy", self.sub_type, s)

    def get_ordered_param_dict(self):
        return Shared.OrderedDict([('genus', self.genus), ('n', self.n), 
                                   ('omegas', self.n_omega), ('degree', self.degree)])


    # since we don't a priori know the number of epsilons, we also need it to compute the partition 
    def get_partition(self, n_epsilon):
        
        # knowing the number of epsilons we can use the previously implemented partition
        V = WOHairyAggregatedGVS(n_components=0,
                                         n_vertices=self.n_vertices,
                                         genus=self.genus,
                                         n=self.n,
                                         n_omega=self.n_omega,
                                         n_epsilon=n_epsilon,
                                         n_double_legs=0)
        
        return V.get_partition()


    # ensure repeatable coloring ----
    def plot_graph(self, G):
        
        GG = Graph(G)  

        n_epsilon = self.get_n_epsilon_from_graph(G)
        
        partition = self.get_partition(n_epsilon)

        assert len(partition) >= self.n + 2

        color_dict = {}

        # inner vertices
        color_dict.update({"gray":partition[0]}) 

        # numbered hairs
        cmap = plt.get_cmap('viridis')
        for i in range(1, self.n+1): 
            assert len(partition[i]) == 1
            interpolation_value = (i - 1) / self.n
            color_dict.update({cmap(interpolation_value):partition[i]})

        # omega vertices
        color_dict.update({"red":partition[self.n + 1]}) 

        # epsilon vertices
        if n_epsilon > 0:
            color_dict.update({"magenta":partition[self.n + 2]}) 

        return GG.plot(vertex_colors=color_dict, vertex_labels=True)



    def is_valid(self):
        # possible values
        l = self.genus >= 1 \
        and self.n >= 0 \
        and self.n_omega >= 0 \
        and self.n_vertices >= 0
        
        return l



    def get_generating_graphs(self):
        #print("Aggregated: get_generating_graphs")
        # Generates all simple graphs with specified number of vertices and edges and at least trivalent vertices.
        if not self.is_valid():
            print("self is not valid")
            return []
        
        # due to symmetry: 
        #   at most one eps-eps double-leg
        #   at most one eps attached to each vertex, or other label: omega or numbered
        max_epsilon = 2 + self.n_vertices + self.n + self.n_omega

        # each connected component has at least one epsilon or omega attached to it
        max_components = self.n_omega + max_epsilon

        for n_epsilon in range(max_epsilon + 1):
            for n_components in range(max_components + 1):
                for n_double_legs in range(n_components + 1):

                    V = WOHairyAggregatedGVS(n_components=n_components,
                                         n_vertices=self.n_vertices,
                                         genus=self.genus,
                                         n=self.n,
                                         n_omega=self.n_omega,
                                         n_epsilon=n_epsilon,
                                         n_double_legs=n_double_legs) 

                    assert V.cohom_degree() == self.degree

                    if V.is_valid():
                        
                        V.build_basis(ignore_existing_files=False)

                        for G in V.get_basis():
                            yield G


    # modified build_basis from GraphVectorSpace class, since we also need to pass n_epsilon to the partition.
    def build_basis(self, progress_bar=False, ignore_existing_files=False, **kwargs):
        """Build the basis of the vector space.

        Create the basis file if the vector space is valid, otherwise skip building a basis. If there exists already
        a basis file rebuild the basis if ignore_existing_file is True, otherwise skip building a basis.
        The basis file contains a list of graph6 strings for canonically labeled graphs building a basis of the vector
        space. The canonical labeling respects the partition of the vertices.

        :param progress_bar: Option to show a progress bar (Default: False).
        :type progress_bar: bool
        :param ignore_existing_files: Option to ignore existing basis file. Ignore existing file and
                rebuild the basis if True, otherwise skip rebuilding the basis file if there exists a basis file already
                (Default: False).
        :type ignore_existing_files: bool
        :param kwargs: Accepting further keyword arguments, which have no influence.
        """
        # print("build basis ", str(self))
        if not self.is_valid():
            # Skip building a basis file if the vector space is not valid.
            return
        if (not ignore_existing_files) and self.exists_basis_file():
            # Skip building a basis file if there exists already one and ignore_existing_file is False.
            return

        generating_list = self.get_generating_graphs()

        desc = 'Build basis: ' + str(self.get_ordered_param_dict())
        # if not progress_bar:
        print(desc)
        basis_set = set()
        # for G in tqdm(generating_list, desc=desc, disable=(not progress_bar)):
        for G in generating_list:

            n_epsilon = self.get_n_epsilon_from_graph(G)
            assert isinstance(n_epsilon, int)

            # For each graph G in the generating list, add the canonical labeled graph6 representation to the basis set
            # if the graph G doesn't have odd automormphisms.
            if self.get_partition(n_epsilon) is None:
                autom_list = G.automorphism_group().gens()
                canonG = G.canonical_label(
                    algorithm=Parameters.canonical_label_algorithm)
            else:
                # The canonical labelling respects the partition of the vertices.
                autom_list = G.automorphism_group(
                    partition=self.get_partition(n_epsilon)).gens()
                canonG = G.canonical_label(partition=self.get_partition(n_epsilon), algorithm=Parameters.canonical_label_algorithm)

            canon6 = canonG.graph6_string()

            if canon6 not in basis_set:
                if not self._has_odd_automorphisms(G, autom_list):
                    basis_set.add(canon6)

        L = list(basis_set)
        L.sort()
        self._store_basis_g6(L)


    # modified graph_to_canon_g6 from GraphVectorSpace class, since we also need to pass n_epsilon to the partition.
    def graph_to_canon_g6(self, graph):
        """Return the graph6 string of the canonically labeled graph and the corresponding permutation sign.

        Labels the sage Graph graph canonically using the sage method for canonical labelling and respecting the
        partition of the vertices.

        :param graph: Graph to be canonically labeled.
        :type graph: Graph
        :return: Tuple containing the graph6 string of the canonically labeled graph and the
            corresponding permutation sign.
        :rtype: tuple(str, int)
        """
        # print("graph_to_canon_g6", graph.graph6_string(), self.get_partition())
        # graph = copy(graph)
        # graph = Graph(graph.graph6_string())

        n_epsilon = self.get_n_epsilon_from_graph(graph)
        assert isinstance(n_epsilon, int)

        canonG, perm_dict = graph.canonical_label(
            partition=self.get_partition(n_epsilon), certificate=True,
            algorithm=Parameters.canonical_label_algorithm)
        
        sgn = self.perm_sign(graph, [perm_dict[j]
                             for j in range(graph.order())])
        assert sgn in {1, -1}

        return (canonG.graph6_string(), sgn)


    def get_work_estimate(self):
        return 0

    @staticmethod
    def compute_deg_min_max(genus, n, n_omega):
        """compute degree lower and upper bounds where the vector space is non-trivial."""

        # degree lower-bound
        # n_vertices = degree - 22 + n_omega - genus + 1 >= 0
        # -> degree >= 22 - n_omega + genus - 1
        deg_min = 22 - n_omega + genus - 1

        # degree upper bound: internal vertices have valcence 3
        # 3 * n_vertices <= 2 * n_edges + (n_hairs - 2*n_double_legs))
        # n_edges = genus + n_vertices + n + n_double_legs - n_hairs - 1
        # -> 3 * n_vertices <= 2 * (genus + n_vertices + n + n_double_legs - n_hairs - 1) + (n_hairs - 2*n_double_legs)
        # -> n_vertices <= 2*genus + 2*n - n_hairs - 2 = 2*genus + n - n_omega - n_epsilon - 2 <= 2*genus + n - n_omega - 2
        # n_vertices = degree - 22 + n_omega - genus + 1
        # -> degree = n_vertices + 22 - n_omega + genus - 1
        # -> degree <= 3*genus + 2*n - n_hairs + 19 - n_omega  
        # -> degree <= 3*genus + 2*n + 19 - 2*n_omega - n - n_epsilon
        # -> degree <= 3*genus + n + 19 - 2*n_omega
        deg_max = 3*genus + n + 19 - 2*n_omega

        excess = 3*(genus - 1) + 2*n - 2*n_omega
        if excess >= 0:
            assert deg_max >= deg_min, "deg_min & deg_max cannot be correct!"

        return deg_min, deg_max


    @staticmethod
    def compute_euler_char(genus, n, n_omega=11):

        euler_char = 0
        excess = 3*(genus - 1) + 2*n - 2*n_omega

        # we add up the euler-characterisitcs of the graph-complex for omega=11,12,13,14,...
        # with each additional omega, the excess decreases by 2
        # further there exist no graphs with negative excess
        # hence we can stop the loop if the excess is negative
        while excess >= 0:

            print("---")
            print("n_omega:", n_omega)
            print("excess:", excess)

            deg_min, deg_max = WOHairyGVS.compute_deg_min_max(genus, n, n_omega)

            euler_char_omega = 0

            for degree in range(deg_min, deg_max + 1):
                
                #print(genus, n, n_omega, degree)
                V = WOHairyGVS(genus=genus, n=n, n_omega=n_omega, degree=degree)

                V.build_basis(ignore_existing_files=False)

                print("degree:", degree, " dimension:", V.get_dimension())
                
                euler_char_omega += (-1)**degree * V.get_dimension()

                if degree > deg_max: assert V.get_dimension() == 0, "deg_max is cannot be correct!"

            print("contribution:", euler_char_omega)
            euler_char += euler_char_omega


            # next iteration
            n_omega += 1
            excess = 3*(genus - 1) + 2*n - 2*n_omega


        return euler_char
    
    
    

        


    



graph_type = "wohairy"


class WOHairyGraphSumVS(GraphVectorSpace.SumVectorSpace):
    """Direct sum of graph vector spaces."""

    def __init__(self, genus_range, n_range, n_omega_range, degree_range):

        self.genus_range = genus_range
        self.n_range = n_range
        self.n_omega_range = n_omega_range
        self.degree_range = degree_range

        vs_list = [WOHairyGVS(genus, n, n_omega, degree) for
                   (genus, n, n_omega, degree) in itertools.product(self.genus_range, self.n_range, self.n_omega_range, self.degree_range)]
        super(WOHairyGraphSumVS, self).__init__(vs_list)

    def get_type(self):
        return 'wohairy graphs'

    def get_ordered_param_range_dict(self):
        return Shared.OrderedDict([('genus', self.genus_range), ('n', self.n_range), ('n_omega', self.n_omega_range), ('degree', self.degree_range)])

    def get_info_plot_path(self):
        s = "info_vector_space_%s" % graph_type
        return os.path.join(Parameters.plots_dir, graph_type, self.sub_type, s)





# ------- Operators --------

class EpsToOmegaGO(SymmetricGraphComplex.SymmetricGraphOperator):
    """Operator that makes one eps into an omega hair."""

    def __init__(self, domain, target):
        super(EpsToOmegaGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return (domain.genus == target.genus 
                and domain.n == target.n
                and domain.n_omega + 1 == target.n_omega 
                and domain.degree - 1 == target.degree)

    @classmethod
    def generate_operator(cls, genus, n, n_omega, degree):
        domain = WOHairyGVS(genus, n, n_omega, degree)
        target = WOHairyGVS(genus, n, n_omega + 1, degree - 1)
        return cls(domain, target)

    def get_matrix_file_path(self):
        s = "epstowD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "epstowD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_type(self):
        return 'eps to omega'


    def operate_on(self, G):
        """
        This is the "dual" operator to the delta_omega operator defined in https://arxiv.org/pdf/2407.16372 on page 4 
        which changes one omega vertex to an epsilon vertex.

        it operates on the graph G by:
        - for each epsilon vertex: replace it with an omega vertex and append the new graph to the image list
        - return the image list

        where replacing the epsilon with an omega vertex works as follows:
        - choose an epsilon vertex
        - note that the first epsilon vertex becomes an omega vertex since target.n_omega = domain.n_omega + 1
        - hence we simply swap the label of the chosen epsilon vertex with the label of the first epsilon vertex
        - the sign is given by the sign of the edge-permutation induced by the relabelling of the vertices
        - note that we do not need to additionaly consider the permutation of the omega-labels, since these are left in place by construction!

        EXAMPLE: "2x omega & 5x epsilon attached to a single vertex"
        notation: "i:0":
        - "i" denotes the "physical" vertex in the graph
        - "0" denotes the corresponding vertex-label

        1.  G1.vertices():  i:0, o1:1, o2:2, e1:3, e2:4, e3:5, e4:6, e5:7
        2.  eps_index = 5 -> eps_vertex e3
        3.  we now want to swap vertex e3 with e1:
            relabeling_perm: [0, 1, 2, 5, 4, 3, 6, 7]
            G2.vertices():  i:0, o1:1, o2:2, e1:5, e2:4, e3:3, e4:6, e5:7

        4.  Now, since each hair-vertex (o1, o2, e1, e2, e3, e4, e5) corresponds to exactly one edge (the hair),
            the induced edge permutation is simply a swap of (i-e3) with (i-e1)
            -> sign = -1
        """

        n_epsilon = len(G.vertices()) - self.domain.n_vertices - self.domain.n - self.domain.n_omega
        assert n_epsilon >= 0

        G1 = copy(G)
        sgn = 1
        image = []

        # label all edges to determine sign later
        Shared.enumerate_edges(G1)

        for i in range(n_epsilon):
            
            omega_eps_cutoff = self.domain.n_vertices + self.domain.n + self.domain.n_omega
            eps_index = omega_eps_cutoff + i
            
            # permute chosen epsilon vertex to the position where it becomes an omega vertex due to the new partition in self.domain
            # if i == 0: then omega_eps_cutoff == eps_index
            G2 = copy(G1)
            if i > 0:
                # relabeling_perm swaps omega_eps_cutoff and eps_index
                relabeling_perm = list(range(omega_eps_cutoff)) + [eps_index] + list(range(omega_eps_cutoff+1, eps_index)) + [omega_eps_cutoff] + list(range(eps_index+1, G1.order()))
                #print(relabeling_perm)
                assert set(relabeling_perm) == set(G2.vertices())
                G2.relabel(relabeling_perm)

            # sign equals the perm_sign of the edge-permutation induced by the relabeling of the vertices
            sgn = Shared.shifted_edge_perm_sign2(G2)
            #print("sgn: ", sgn)
            image.append((G2, sgn))

        return image


    def get_work_estimate(self):
        return 0



class EpsToOmegaD(GraphOperator.Differential):

    def __init__(self, sum_vector_space):
        """Initialize the eps to omega differential with the underlying sum vector space."""

        super(EpsToOmegaD, self).__init__(sum_vector_space,
                                          EpsToOmegaGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'EpsToOmega'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_epstoomega_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)
    







class ContractEdgesGO(SymmetricGraphComplex.SymmetricGraphOperator):
    """Contract edges graph operator."""

    def __init__(self, domain, target):
        super(ContractEdgesGO, self).__init__(domain, target)

    @staticmethod
    def is_match(domain, target):
        return (domain.genus == target.genus 
                and domain.n == target.n
                and domain.n_omega == target.n_omega 
                and domain.degree - 1 == target.degree)

    @classmethod
    def generate_operator(cls, genus, n, n_omega, degree):
        domain = WOHairyGVS(genus, n, n_omega, degree)
        target = WOHairyGVS(genus, n, n_omega, degree-1)
        return cls(domain, target)


    def get_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_rank_file_path(self):
        s = "contractD%d_%d_%d_%d_rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.data_dir, graph_type, s)

    def get_ref_matrix_file_path(self):
        s = "contractD%d_%d_%d_%d.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_ref_rank_file_path(self):
        s = "contractD%d_%d_%d_%d.txt.rank.txt" % self.domain.get_ordered_param_dict().get_value_tuple()
        return os.path.join(Parameters.ref_data_dir, graph_type, s)

    def get_type(self):
        return 'contract edges'


    def operate_on(self, G):
        """
        This is the "dual" operator to the splitting-vertices operator delta_s defined in https://arxiv.org/pdf/2407.16372 on page 4.

        It Operates on the graph G by going though all edges:
        For each edge it first checks whether it can be contracted, i.e. whether it connects at least one internal vertex and is not connected to a numbered hair-vertex.
        Then it contracts the edge and unifies the adjacent vertices.
        """
        
        n_vertices = self.domain.n_vertices
        n = self.domain.n
        n_omega = self.domain.n_omega
        n_epsilon = self.domain.get_n_epsilon_from_graph(G)

        image = []
        for (i, e) in enumerate(G.edges(labels=False, sort=True)):

            #print("contracting edge (u,v) =", e, "###########")
            
            (u, v) = e
            sgn = (-1)**i

            # ensure u<v (this should be always true anyway actually)
            assert u < v

            # only edges connected to at least one internal vertex, and not connected to a numbered hair-vertex can be contracted

            # both u and v are not internal vertices
            if u >= n_vertices: continue

            # u or v are numbered vertices        
            if u >= n_vertices and u < n_vertices + n: continue
            if v >= n_vertices and v < n_vertices + n: continue

            # label all edges to determine sign later
            G1 = copy(G)
            Shared.enumerate_edges(G1)


            # contracting the edge (u,v) ---

            # if both u and v are internal vertices
            if v < n_vertices:
                
                #print("CASE: inner_vertex - inner_vertex")

                assert u < n_vertices

                # contract edge by merging vertices
                G1.merge_vertices([v, u])

                # if we have too few edges, some double edges have been created => zero
                if G.size() - G1.size() > 1: continue

                # relabel vertices back to 0,1,2,...,k
                G1.relabel(range(len(G1.vertices())), inplace=True)
                
                # find edge permutation sign
                sgn *= Shared.shifted_edge_perm_sign2(G1)
                # print("sgn3_",sgn)
                image.append((G1, sgn))
                # image.append((Graph(G1.graph6_string()), sgn))
                # print("hmm0:", G.graph6_string(), G1.graph6_string())


            # if u is an internal vertex and v is an epsilon-vertex
            elif u < n_vertices \
                and v >= n_vertices + n + n_omega:
                
                #print("CASE: inner_vertex - epsilon")

                G1.delete_vertex(v)

                # split u into separate epsilon-vertices
                for w in G1.neighbors(u):

                    #print("picking neighbour w =", w, " ---")
                    
                    u_w_label = G1.edge_label(u, w)
                    G1.delete_edge(u, w)

                    new_eps_label = max(G1.vertices()) + 1
                    assert not new_eps_label in G1.vertices()

                    G1.add_vertex(new_eps_label)
                    G1.add_edge(w, new_eps_label, u_w_label)

                G1.delete_vertex(u) 

                # relabel vertices back to 0,1,2,...,k
                G1.relabel(range(len(G1.vertices())), inplace=True)

                sgn *= Shared.shifted_edge_perm_sign2(G1)
                image.append((G1, sgn))


            # if u is an internal vertex and v is a omega-vertex
            # the second vertex is now an omega-vertex, so we need to merge the vertex with the eps vertex
            # after reconnecting one of the edges to omega
            # we assume that u != eps, because eps-omega-edges cannot be contracted
            elif u < n_vertices \
                and v >= n_vertices + n \
                and v < n_vertices + n + n_omega:
                
                #print("CASE: inner_vertex - omega")

                G1.delete_vertex(v)

                # pick vertex w which will be connected to omega
                for (j, w) in enumerate(G1.neighbors(u)):

                    #print("picking neighbour w =", w, " ---")
                    G2 = copy(G1)
                    sgn2 = sgn

                    u_w_label = G2.edge_label(u, w)
                    G2.delete_edge(u, w)

                    # why v is convenient:
                    # - v was the label of some omega-vertex 
                    #   ->  it will again correspond to an omega vertex in the new partition in self.target
                    #       as long as we insure to add the new epsilons after it
                    # - v has been deleted and is hence a "free" vertex-label
                    new_omega_label = v
                    assert not new_omega_label in G2.vertices()

                    G2.add_vertex(new_omega_label)
                    G2.add_edge(w, new_omega_label, u_w_label)

                    # all other vertices s will be connected to epsilons
                    n_new_eps = len(G2.neighbors(u))
                    for s in G2.neighbors(u): 
                        
                        #print("picking neighbour s =", s, " ---")
                    
                        u_s_label = G2.edge_label(u, s)
                        G2.delete_edge(u, s)

                        new_eps_label = max(G2.vertices()) + 1
                        assert not new_eps_label in G2.vertices()

                        G2.add_vertex(new_eps_label)
                        G2.add_edge(s, new_eps_label, u_s_label)
                    
                    G2.delete_vertex(u) 

                    # relabel vertices back to 0,1,2,...,k
                    G2.relabel(range(len(G2.vertices())), inplace=True)
                    
                    sgn2 *= Shared.shifted_edge_perm_sign2(G2)

                    # sanity-checks
                    n_eps_in_target = self.target.get_n_epsilon_from_graph(G2)
                    assert n_eps_in_target == n_epsilon + n_new_eps
                    assert G2.order() == self.target.n_vertices + self.target.n + self.target.n_omega + n_eps_in_target

                    image.append((G2, sgn2))


        return image
    
    
    def get_work_estimate(self):
        return 0


class ContractEdgesD(GraphOperator.Differential):
    """Contract edges differential."""

    def __init__(self, sum_vector_space):
        """Initialize the contract edges differential with the underlying sum vector space."""

        super(ContractEdgesD, self).__init__(sum_vector_space,
                                             ContractEdgesGO.generate_op_matrix_list(sum_vector_space))

    def get_type(self):
        return 'contract edges'

    def get_cohomology_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)

    def get_cohomology_web_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "cohomology_dim_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.web_dir, graph_type, sub_type, s)

    def get_info_plot_path(self):
        sub_type = self.sum_vector_space.sub_type
        s = "info_contract_D_%s" % (graph_type)
        return os.path.join(Parameters.plots_dir, graph_type, sub_type, s)




# ------- Graph Complex --------
class WOHairyGC(GraphComplex.GraphComplex):

    def __init__(self, genus_range, n_range, omega_range, degree_range, differentials):
        """Initialize the graph complex."""
        
        self.genus_range = genus_range
        self.n_range = n_range
        self.omega_range = omega_range
        self.degree_range = degree_range

        sum_vector_space = WOHairyGraphSumVS(
            self.genus_range, self.n_range, self.omega_range, self.degree_range)
        
        differential_list = []

        if not set(differentials).issubset(['contract', 'contract_iso', 'epstoomega', 'epstoomega_iso']):
            raise ValueError(
                "Differentials for hairy graph complex: 'contract'")
        
        contract_edges_dif = ContractEdgesD(sum_vector_space)
        epstoomega_dif = EpsToOmegaD(sum_vector_space)

        if 'contract' in differentials:
            differential_list.append(contract_edges_dif)
        
        if 'epstoomega' in differentials:
            differential_list.append(epstoomega_dif)
        super(WOHairyGC, self).__init__(sum_vector_space, differential_list)


    def __str__(self):
        return '<%s graph complex with %s>' % (graph_type, str(self.sub_type))
    


    @staticmethod
    def DSquareTest_single(operator, genus, n, n_omega=11):
        """Testing the D^2=0 property of either ContractEdges or EpsToOmega"""

        assert operator in ['contract', 'epstoomega']

        excess = 3*(genus - 1) + 2*n - 2*n_omega
        assert excess >= 0
        n_omega_max = n_omega + (excess // 2)

        deg_min, deg_max = WOHairyGVS.compute_deg_min_max(genus, n, n_omega)

        GC = WOHairyGC(genus_range=[genus], 
                n_range=[n], 
                omega_range=range(n_omega, n_omega_max + 1), 
                degree_range=range(deg_min, deg_max + 1), 
                differentials=[operator])

        GC.build_basis(ignore_existing_files=False)
        GC.build_matrix(ignore_existing_files=True)


        for dif in GC.operator_collection_list:

            succ_l = 0
            triv_l = 0
            
            for (op1, op2) in itertools.permutations(dif.op_matrix_list, 2):

                if op1.get_target() == op2.get_domain():
                    # A composable pair is found

                    pair = (op1, op2)
                    res = dif._square_zero_test_for_pair(pair, eps=Parameters.square_zero_test_eps)

                    if res == 'triv':
                        triv_l += 1
                    elif res == 'succ':
                        succ_l += 1
                    else:
                        assert False, "d_square zero test failed for %s" % str(pair)

            print("trivial success:", triv_l)
            print("success:", succ_l)


    @staticmethod
    def Anticomm_Test_single(genus, n, n_omega=11, eps=Parameters.square_zero_test_eps):
        """Testing the anti-commutativity of ContractEdges and EpsToOmega operators"""

        excess = 3*(genus - 1) + 2*n - 2*n_omega
        assert excess >= 0
        n_omega_max = n_omega + (excess // 2)

        deg_min, deg_max = WOHairyGVS.compute_deg_min_max(genus, n, n_omega)

        GC = WOHairyGC(genus_range=[genus], 
                n_range=[n], 
                omega_range=range(n_omega, n_omega_max + 1), 
                degree_range=range(deg_min, deg_max + 1), 
                differentials=['contract', 'epstoomega'])

        GC.build_basis(ignore_existing_files=False)
        GC.build_matrix(ignore_existing_files=False)

        print(GC.operator_collection_list)
        assert len(GC.operator_collection_list) == 2

        dif_1 = GC.operator_collection_list[0]
        dif_2 = GC.operator_collection_list[1]
        op1_matrix_list = dif_1.op_matrix_list
        op2_matrix_list = dif_2.op_matrix_list
        
        for (op1_1, op1_2, op2_1, op2_2) in itertools.product(op1_matrix_list, op1_matrix_list, op2_matrix_list, op2_matrix_list):
            if op1_1.get_domain() == op2_2.get_domain() and op1_2.get_target() == op2_1.get_target() \
                and op1_1.get_target() == op2_1.get_domain() and op2_2.get_target() == op1_2.get_domain():
                if (not (op1_1.is_valid() and op1_2.is_valid() and op2_1.is_valid() and op2_2.is_valid())) \
                    or ((op1_1.is_trivial() or op2_1.is_trivial()) and (op1_2.is_trivial() or op2_2.is_trivial())):
                    #print('trivial success')
                    pass

                else: 
                    M1_1 = op1_1.get_matrix()
                    M1_2 = op1_2.get_matrix()
                    M2_1 = op2_1.get_matrix()
                    M2_2 = op2_2.get_matrix()

                    product_1 = M2_1 * M1_1
                    product_2 = M1_2 * M2_2

                    if Shared.matrix_norm(product_1 - product_2) < eps:
                        print('success')
                    else:
                        assert False, 'Anitcomm-Test failed!'


    @staticmethod
    def compute_cohomology_dim(degree, genus, n, n_omega=11, prev_r2=None):
        """compute the cohomology dimension for a given (genus, n) pair"""
        
        if prev_r2 != None: 
            assert isinstance(prev_r2, int)
            assert prev_r2 >= 0


        D_list = []

        # the r2 from the previous degree is equal to the r1 of the current degree
        # hence, if we already have it, we can skip the computation of r1
        if prev_r2 == None: 
            D_Cont_deg = ContractEdgesGO.generate_operator(degree=degree, genus=genus, n=n, n_omega=n_omega)
            D_eps_deg = EpsToOmegaGO.generate_operator(degree=degree, genus=genus, n=n, n_omega=n_omega)
            D_list.append(D_Cont_deg)
            D_list.append(D_eps_deg)

        D_Cont_deg_p1 = ContractEdgesGO.generate_operator(degree=degree+1, genus=genus, n=n, n_omega=n_omega)
        D_eps_deg_p1 = EpsToOmegaGO.generate_operator(degree=degree+1, genus=genus, n=n, n_omega=n_omega)
        D_list.append(D_Cont_deg_p1)
        D_list.append(D_eps_deg_p1)


        for D in D_list:
            D.domain.build_basis(skip_existing_files=False)
            D.target.build_basis(skip_existing_files=False)

            D.build_matrix(skip_existing_files=False)


        if prev_r2 == None: 
            D_Cont_deg_mat = D_Cont_deg.get_matrix()
            D_eps_deg_mat = D_eps_deg.get_matrix()

        D_Cont_deg_p1_mat = D_Cont_deg_p1.get_matrix()
        D_eps_deg_p1_mat = D_eps_deg_p1.get_matrix()

        if prev_r2 == None: 
            assert D_Cont_deg.domain == D_eps_deg.domain
            deg_double_mat = D_Cont_deg_mat.stack(D_eps_deg_mat)

        assert D_Cont_deg_p1.domain == D_eps_deg_p1.domain
        deg_p1_double_mat = D_Cont_deg_p1_mat.stack(D_eps_deg_p1_mat) 
        

        d = WOHairyGVS(genus=genus, n=n, n_omega=n_omega, degree=degree).get_dimension()

        print("computing ranks for: degree =", degree, ", dimension =", d)

        if prev_r2 == None: 
            r1 = deg_double_mat.rank()
        else:
            print("r1 = r2 from previous degree")
            r1 = prev_r2
        
        print("computing r2")
        r2 = deg_p1_double_mat.rank()
        print("computing r3")
        r3 = D_eps_deg_p1_mat.rank()

        cohom_dim = d - r1 - r2 + r3


        assert cohom_dim >= 0
        if cohom_dim > 0: print("k =", degree, ":", cohom_dim)

        return cohom_dim, r2


    @staticmethod
    def max_basis_dimension(genus, n, n_omega=11):
        """
        used for training the model in WOHairyGC_get_dimension_predictor.py
        which again lets us obtain: max_basis_dimension_estimate()
        """
 
        excess = 3*(genus - 1) + 2*n - 2*n_omega
        if excess < 0: return 0

        deg_min, deg_max = WOHairyGVS.compute_deg_min_max(genus, n, n_omega)

        max_dim = 0

        for degree in range(deg_min, deg_max + 1):
            
            #print(genus, n, n_omega, degree)
            V = WOHairyGVS(genus=genus, n=n, n_omega=n_omega, degree=degree)

            V.build_basis(ignore_existing_files=False)

            max_dim = max(max_dim, V.get_dimension())
            
            if degree > deg_max: assert V.get_dimension() == 0, "deg_max is cannot be correct!"

        return max_dim


    @staticmethod
    # model obtained from WOHairyGC_get_dimension_predictor.py
    def max_basis_dimension_estimate(genus, n):
        n_exp = 1.7
        coefficients = [2.84547639, 0.33797069, 0.24190675]
        intercept = -24.538197700729306
        return int(np.exp(coefficients[0]*genus + coefficients[1]*(n**n_exp) + coefficients[2]*genus*n + intercept))


    @staticmethod
    def compute_cohomology_dim_all(g_max=20, n_max=20, n_omega=11):
        """
        computes the cohomology dimensions for all (genus, n) pairs by work-estimate
        saves the results in a the csv file 'WOHairy_CohomologyDimensions.csv'
        """

        jobs = []
        for genus in range(1, g_max + 1):
            for n in range(n_max + 1):
                basis_estimate = WOHairyGC.max_basis_dimension_estimate(genus, n)
                if basis_estimate < 1000000:
                    jobs.append((genus, n, basis_estimate))

        table = [["?" for _ in range(n_max+2)] for _ in range(1,g_max+2)]
        table[0][0] = "g / n"
        for n in range(0, n_max+1): 
            table[0][n+1] = str(n)
        for g in range(1, g_max+1):
            table[g][0] = str(g)

        # sort jobs by basis estimate
        jobs.sort(key=lambda x: x[2])
        print("g_n's sorted by basis estimate")
        print(jobs)

        # try to load the table from a file if it already exists
        # this prevents recomputation of ranks
        try:
            with open("WOHairy_CohomologyDimensions.csv", "r", newline="") as file:
                reader = csv.reader(file)
                table = [row for row in reader]
        except FileNotFoundError:
            pass  # no table file yet

        already_computed_g_n = []
        for g in range(1,g_max+1):
            for np1 in range(1,n_max+2):
                if table[g][np1] != "?":
                    already_computed_g_n.append((g, np1-1))
        print("already_computed_g_n")
        print(already_computed_g_n)


        for genus, n, basis_estimate in jobs:

            if (genus, n) in already_computed_g_n:
                print("")
                print("Already computed (genus, n) = ", (genus, n))
                print("max_basis_dim:", WOHairyGC.max_basis_dimension(genus, n), " / ", WOHairyGC.max_basis_dimension_estimate(genus, n), "(estimate)")
                continue
            
            print("")
            print("(genus, n) = ", (genus, n))
            excess = 3*(genus - 1) + 2*n - 2*n_omega
            print("excess:", excess)

            print("max_basis_dim_estimate:", WOHairyGC.max_basis_dimension_estimate(genus, n))

            deg_min, deg_max = WOHairyGVS.compute_deg_min_max(genus, n, n_omega)

            degree_range = range(deg_min, deg_max+3)

            non_zero_dim_list = []

            r2 = None
            for degree in degree_range:

                cohom_dim, r2 = WOHairyGC.compute_cohomology_dim(degree=degree, genus=genus, n=n, n_omega=n_omega, prev_r2=r2)

                if cohom_dim > 0:
                    non_zero_dim_list.append((degree, cohom_dim))

            if len(non_zero_dim_list) > 0:
                table[genus][n+1] = ""
                for i, (degree, cohom_dim) in enumerate(non_zero_dim_list):
                    table[genus][n+1] += "k=" + str(degree) + ": " + str(cohom_dim)

                    if i < len(non_zero_dim_list) - 1:
                        table[genus][n+1] += "\n"
            else:
                table[genus][n+1] = "0"

            print("max_basis_dim:", WOHairyGC.max_basis_dimension(genus, n), " / ", WOHairyGC.max_basis_dimension_estimate(genus, n), "(estimate)")
            print("Starting to write 'WOHairy_CohomologyDimensions.csv'")

            # Save the table to a CSV file
            with open("WOHairy_CohomologyDimensions.csv", "w", newline="") as file:
                writer = csv.writer(file)
                writer.writerows(table)

            print("Table saved to 'WOHairy_CohomologyDimensions.csv'")








    





    



