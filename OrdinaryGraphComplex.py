import os
from sage.all import *
import GraphVectorSpace as GVS

reload(GVS)

class OrdinaryGraphVectorSpace(GVS.GraphVectorSpace):
    dataDir = "./data"
    dataDirOdd = dataDir + "/ordinary/oddedge/"
    dataDirEven = dataDir + "/ordinary/evenedge/"
    imgBaseDir = "img/"

    def __init__(self, nVertices, nLoops, evenEdges=True):
        self.nVertices = nVertices
        self.nLoops = nLoops
        self.evenEdges = evenEdges

    def file_name(self):
        dataDir = OrdinaryGraphVectorSpace.dataDirEven if self.evenEdges else OrdinaryGraphVectorSpace.dataDirOdd
        s = "gra%d_%d.g6" % (self.nVertices, self.nLoops)
        return os.path.join(dataDir, s)

    def svg_dir(self):
        dataDir = OrdinaryGraphVectorSpace.dataDirEven if self.evenEdges else OrdinaryGraphVectorSpace.dataDirOdd
        s = "imgs%d_%d/" % (self.nVertices, self.nLoops)
        return os.path.join(dataDir, OrdinaryGraphVectorSpace.imgBaseDir, s)

    def color_counts(self):
        return None # no coloring

    def valid(self):
        nEdges = self.nLoops + self.nVertices -1
        # at least trivalent, and simple
        return  (3*self.nVertices <= 2*nEdges) and self.nVertices > 0 and self.nLoops >= 0 and nEdges <= self.nVertices*(self.nVertices-1)/2

    def _generating_graphs(self):
        # generate List of unmarked graphs
        graphList = self._list_graphs(self.nVertices, self.nLoops)
        #graphList = slef._listGraphs(self.nVertices, self.nLoops, false)
        return graphList

    def _list_graphs(self, nVertices, nLoops, onlyonevi=True):
        """
        creates a list of simple 1vi graphs with at least trivalent vertices
        """
        nEdges = nLoops + nVertices - 1
        if (3 * nVertices > 2 * nEdges) or (nEdges > nVertices * (nVertices - 1) / 2):
            # impossible
            return []
        gL = list(graphs.nauty_geng(("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (nVertices, nEdges, nEdges)))
        return gL

    def _perm_sign(self, G, p):
        nVert, nLoops, evenEdges = (self.nVertices, self.nLoops, self.evenEdges)
        nEdges = nLoops + nVert - 1
        if evenEdges:
            # The sign is (induced sign on vertices) * (induced sign edge orientations)
            pp=[j+1 for j in p]
            sgn = Permutation(pp).signature()
            for (u, v, ignored) in G.edges():
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
                u,v = e
                G1.set_edge_label(u,v,j)

            # we permute the graph, and read of the new labels
            G1.relabel(perm=p,inplace=True)
            pp = [j+1 for u,v,j in G1.edges()]
            return Permutation(pp).signature()

    def canonical(self, graph, colorData=None):
        return graph.canonical_label()

    def work_estimate(self):
        # give estimate of number of graphs
        nEdges = self.nLoops + self.nVertices - 1
        n = self.nVertices
        return binomial((n*(n-1))/2, nEdges) / factorial(n)

