import os
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator as GO

reload(GVS)

dataDir = "./data"
dataDirOdd = dataDir + "/ordinary/oddedge/"
dataDirEven = dataDir + "/ordinary/evenedge/"
imgBaseDir = "img/"

class OrdinaryGraphVectorSpace(GVS.GraphVectorSpace):

    def __init__(self, nVertices, nLoops, evenEdges=True):
        self.nVertices = nVertices
        self.nLoops = nLoops
        self.evenEdges = evenEdges

    def file_name(self):
        directory = dataDirEven if self.evenEdges else dataDirOdd
        s = "gra%d_%d.g6" % (self.nVertices, self.nLoops)
        return os.path.join(directory, s)

    def svg_dir(self):
        directory = dataDirEven if self.evenEdges else dataDirOdd
        s = "imgs%d_%d/" % (self.nVertices, self.nLoops)
        return os.path.join(directory, imgBaseDir, s)

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

    def work_estimate(self):
        # give estimate of number of graphs
        nEdges = self.nLoops + self.nVertices - 1
        n = self.nVertices
        return binomial((n*(n-1))/2, nEdges) / factorial(n)


# -----  Contraction operator --------

class ContractDOrdinary(GO.GraphOperator):

    def __init__(self, nVertices, nLoops, evenEdges=True):
        # source vector space:
        self.nVertices = nVertices
        self.nLoops = nLoops
        self.evenEdges = evenEdges

    def file_name(self):
        directory = directory = dataDirEven if self.evenEdges else dataDirOdd
        s = "contractD%d_%d.txt" % (self.nVertices, self.nLoops)
        return os.path.join(directory, s)

    def domain(self):
        return OrdinaryGraphVectorSpace(self.nVertices, self.nLoops, self.evenEdges)

    def target(self):
        return OrdinaryGraphVectorSpace(self.nVertices-1, self.nLoops, self.evenEdges)

    def _operate_on(self,graph):

        image=[]
        for (u,v) = edges(G)
        # permute u,v to position 1 and 2
        p = collect(1:self.nVertices)
        p[1] = u
        p[2] = v
        idx = 3
        for j = 1:self.nVertices
            if j == u || j== v
                continue
            else
                p[idx] = j
                idx +=1

        pp = invPermutation(p)
        #println(pp)
        sgn = get_perm_sign(vs,G, pp)
        GG = permuteGraph(G,pp)

        # now delete the first edge
        remove_edge!(GG,1,2) #.... done later
        # ... and join 0 and 1 and fix labelings
        p = vcat([1],collect(1:self.nVertices-1))
        if !self.evenEdges
          sgn *= get_perm_sign(vs,GG, p)  # TODO: no good to call get_perm_sign with non-permutation if evenEdges
        GG = permuteGraph(GG,p)

        # finally delete the last vertex
        # remove_vertex(GG, nVert) ... by creating a new graph
        GGG = small_graph(self.nVertices-1)
        for (uu,vv) = edges(GG)
          if uu != vv
            add_edge!(GGG, uu, vv)

        push!(ret, (GGG, sgn))
        return ret

    def work_estimate(self):
        # give estimate of number of graphs
        domain = self.domain()
        nEdges = domain.nLoops + domain.nVertices -1
        return domain.work_estimate * nEdges