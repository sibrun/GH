import os
from sage.all import *
from GraphVectorSpace import GraphVectorSpace


class OrdinaryGraphVectorSpace(GraphVectorSpace):

    ordinaryDataDirOdd = dataDir + "/ordinarydata/oddedge/"
    ordinaryDataDirEven = dataDir + "/ordinarydata/evenedge/"
    imgBaseDir = "img/"

    def __init__(self, nVertices=0, nLoops=0, evenEdges=True):
        self.nVertices = nVertices
        self.nLoops = nLoops
        self.evenEdges = evenEdges

    def get_file_name(self):
        dataDir = OrdinaryGraphVectorSpace.ordinaryDataDirEven if self.evenEdges else OrdinaryGraphVectorSpace.ordinaryDataDirOdd
        s = "gra" + self.nVertices + "_" + self.nLoops + ".g6"
        return os.path.join(dataDir, s)

    def get_svg_dir(self):
        dataDir = OrdinaryGraphVectorSpace.ordinaryDataDirEven if self.evenEdges else OrdinaryGraphVectorSpace.ordinaryDataDirOdd
        s = "imgs%d_%d/" % (self.nVertices, self.nLoops)
        return os.path.join(dataDir, OrdinaryGraphVectorSpace.imgBaseDir, s)

    def get_color_counts(self):
        return None # no coloring

    def is_valid(self):
        nEdges = self.nLoops + self.nVertices -1
        # at least trivalent, and simple
        return  (3*self.nVertices <= 2*nEdges) and self.nVertices > 0 and self.nLoops >= 0 and nEdges <= self.nVertices*(self.nVertices-1)/2

    def get_generating_graphs(self):
        # generate List of unmarked graphs
        #println(self.nVertices)
        graphList = self._listGraphs(self.nVertices, self.nLoops)
        #LL = listG(self.nVertices, self.nLoops, false)
        #LL=[]
        return graphList

    def _listGraphs(nVertices, nLoops, onlyonevi=True):
        """
        creates a list of simple 1vi graphs with at least trivalent vertices
        """
        nEdges = nLoops + nVertices - 1
        if (3 * nVertices > 2 * nEdges) || (nEdges > nVertices * (nVertices - 1) / 2):
            # impossible
            return []

        gL = list(graphs.nauty_geng(("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (nVertices, nEdges, nEdges)))
        return gL

    def get_perm_sign(self, G, p):
        nVert, nLoops, evenEdges = (self.nVertices, self.nLoops, self.evenEdges)
        nEdges = nLoops + nVert - 1
        if evenEdges:
            # The sign is (induced sign on vertices) * (induced sign edge orientations)
            sgn = p.signature()
            for (u,v,ignored) in G.edges():
                # we assume the edge is always directed from the larger to smaller index
                if (u < v and p[u] > p[v]) or (u > v and p[u] < p[v]):
                    sgn *= -1

            return sgn
        else:
            # The sign is (induced sign of the edge permutation)
            # we assume the edges are always lex ordered
            # for the computation we use that edges(G) returns the edges in lex ordering
            #pp = range(nEdges)
            ##G2 = permuteGraph(G,p)

            # we first label the edges on a copy of G lexicographically
            G1 = copy(G)
            for (j,e) in enumerate(G1.edges(label=False)):
                u,v = e
                G1.set_edge_label(u,v,j)

            # we permute the graph, and read of the new labels
            G1.relabel(p,inplace=True)
            pp = [j+1 for u,v,j in G1.edges()]

            #println(pp)
            return Permutation(pp).sign()

    def get_work_estimate(self):
        # give estimate of number of graphs
        nEdges = self.nLoops + self.nVertices -1
        n = self.nVertices
        return binomial((n*(n-1))//2, nEdges) / factorial(n)



#"""Converts the graph to a graphviz dot format string.
#   This method is used only for visualization, not for computation."""
#function get_dot(self::OrdinaryGraphVectorSpace, G)
#    return render_to_dot(G)
#end

