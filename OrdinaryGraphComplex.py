import os
import itertools
import operator
import logging
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO
import GraphComplex as GC
import Shared as SH

reload(GVS)
reload(GO)
reload(GC)
reload(SH)

data_dir = "./data"
data_ref_dir = "./data_ref"
type_dir = "ordinary"
sub_dir_odd = "oddedge"
sub_dir_even = "evenedge"
image_directory = "img"


class OrdinaryGVS(GVS.GraphVectorSpace):

    def __init__(self, n_vertices, n_loops, even_edges=True, header_ref=False):
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1
        super(OrdinaryGVS,self).__init__(header_ref=header_ref)

    def _set_file_name(self, ref=False):
        s0 = data_ref_dir if ref else data_dir
        s1 = sub_dir_even if self.even_edges else sub_dir_odd
        s2 = "gra%d_%d.g6" % (self.n_vertices, self.n_loops)
        return os.path.join(s0, type_dir, s1, s2)

    def _set_validity(self):
        return (3 * self.n_vertices <= 2 * self.n_edges) and self.n_vertices > 0 and self.n_loops >= 0 and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2

    def _set_work_estimate(self):
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)

    def params_to_string(self):
        return "n_vertices=%d, n_loops=%d, even_edges=%d" % (self.n_vertices, self.n_loops, self.even_edges)

    def _generating_graphs(self, onlyonevi=True):
        if (3 * self.n_vertices > 2 * self.n_edges) or (self.n_edges > self.n_vertices * (self.n_vertices - 1) / 2):
            # impossible
            return []
        return list(graphs.nauty_geng(("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (self.n_vertices, self.n_edges, self.n_edges)))

    def perm_sign(self, G, p):
        if self.even_edges:
            # The sign is (induced sign on vertices) * (induced sign edge orientations)
            sgn = SH.Perm(p).sign()
            for (u, v) in G.edges(labels=False):
                # we assume the edge is always directed from the larger to smaller index
                if (u < v and p[u] > p[v]) or (u > v and p[u] < p[v]):
                    sgn *= -1
            return sgn
        else:
            # The sign is (induced sign of the edge permutation)
            # we assume the edges are always lex ordered
            # for the computation we use that G.edges() returns the edges in lex ordering
            # we first label the edges on a copy of G lexicographically
            G1 = copy(G)
            for (j,e) in enumerate(G1.edges(labels=False)):
                (u, v) = e
                G1.set_edge_label(u, v, j)

            # we permute the graph, and read of the new labels
            G1.relabel(p, inplace=True)
            return SH.Perm([j for (u, v, j) in G1.edges()]).sign()


# -----  Contraction operator --------
class ContractGO(GO.GraphOperator):

    def __init__(self, n_vertices, n_loops, even_edges=True, header_ref=False, skip_if_no_basis=True):
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges

        domain = OrdinaryGVS(self.n_vertices, self.n_loops, even_edges=self.even_edges)
        target = OrdinaryGVS(self.n_vertices - 1, self.n_loops, even_edges=self.even_edges)

        super(ContractGO, self).__init__(domain, target, header_ref=header_ref, skip_if_no_basis=skip_if_no_basis)

    def _set_file_name(self, ref=False):
        s0 = data_ref_dir if ref else data_dir
        s1 = sub_dir_even if self.even_edges else sub_dir_odd
        s2 = "contractD%d_%d.txt" % (self.n_vertices, self.n_loops)
        return os.path.join(s0, type_dir, s1, s2)

    def _set_work_estimate(self):
        return self.domain.work_estimate * self.domain.n_edges

    def params_to_string(self):
        return "n_vertices=%d, n_loops=%d, even_edges=%d" % (self.n_vertices, self.n_loops, self.even_edges)

    def _operate_on(self,G):
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            r = range(0,self.domain.n_vertices)
            p = list(r)
            p[0] = u
            p[1] = v
            idx = 2
            for j in r:
                if j == u or j== v:
                    continue
                else:
                    p[idx] = j
                    idx +=1

            pp = SH.Perm(p).inverse()
            sgn = self.domain.perm_sign(G, pp)
            G1 = copy(G)
            G1.relabel(pp, inplace=True)

            for (j, ee) in enumerate(G1.edges(labels=False)):
                a, b = ee
                G1.set_edge_label(a,b,j)
            previous_size = G1.size()
            G1.merge_vertices([0,1])
            if (previous_size - G1.size()) != 1:
                continue
            G1.relabel(list(range(0,G1.order())), inplace = True)
            if not self.domain.even_edges:
                p = [j for (a, b, j) in G1.edges()]
                sgn *= Permutation(p).signature()
            image.append((G1, sgn))
        return image


# ----- Ordinary Graph Complex --------
class OrdinaryGC(GC.GraphComplex):
    def __init__(self, v_range, l_range, even_range, header_ref=False, skip_if_no_basis=True, delete_old=False):
        self.v_range = v_range
        self.l_range = l_range
        self.even_range = even_range
        self.header_ref = header_ref
        self.skip_if_no_basis = skip_if_no_basis
        super(OrdinaryGC, self).__init__(delete_old=delete_old)

    def create_vs(self):
        for (v, l, even) in itertools.product(self.v_range, self.l_range, self.even_range):
            self.vs_list.append(OrdinaryGVS(v, l, even_edges=even, header_ref=self.header_ref))

    def create_op(self):
        for (v, l, even) in itertools.product(self.v_range, self.l_range, self.even_range):
            self.op_list.append(ContractGO(v, l, even_edges=even, header_ref=self.header_ref, skip_if_no_basis=self.skip_if_no_basis))

