import os
from sage.all import *
import GraphVectorSpace as GVS
import GraphOperator

reload(GVS)
reload(GraphOperator)

data_directory = "./data"
data_directory_odd = data_directory + "/ordinary/oddedge/"
data_directory_even = data_directory + "/ordinary/evenedge/"
imgBaseDir = "img/"


class OrdinaryGVS(GVS.GraphVectorSpace):

    def __init__(self, n_vertices, n_loops, even_edges=True):
        self.n_vertices = n_vertices
        self.n_loops = n_loops
        self.even_edges = even_edges
        self.n_edges = self.n_loops + self.n_vertices - 1

        directory = data_directory_even if self.even_edges else data_directory_odd
        s = "gra%d_%d.g6" % (self.n_vertices, self.n_loops)
        file_name = os.path.join(directory, s)

        valid = (3 * self.n_vertices <= 2 * self.n_edges) and self.n_vertices > 0 and self.n_loops >= 0 and self.n_edges <= self.n_vertices * (self.n_vertices - 1) / 2

        super(OrdinaryGVS,self).__init__(valid, file_name)

    def _generating_graphs(self):
        # generate List of unmarked graphs
        graphList = self._list_graphs()
        #graphList = slef._listGraphs(self.nVertices, self.nLoops, false)
        return graphList

    def _list_graphs(self, onlyonevi=True):
        """
        creates a list of simple 1vi graphs with at least trivalent vertices
        """
        if (3 * self.n_vertices > 2 * self.n_edges) or (self.n_edges > self.n_vertices * (self.n_vertices - 1) / 2):
            # impossible
            return []
        gL = list(graphs.nauty_geng(("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (self.n_vertices, self.n_edges, self.n_edges)))
        return gL

    def perm_sign(self, G, p):
        if self.even_edges:
            # The sign is (induced sign on vertices) * (induced sign edge orientations)
            sgn = Permutation([j+1 for j in p]).signature()
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
            return Permutation([j+1 for (u, v, j) in G1.edges()]).signature()

    def get_work_estimate(self):
        # give estimate of number of graphs
        return binomial((self.n_vertices * (self.n_vertices - 1)) / 2, self.n_edges) / factorial(self.n_vertices)


# -----  Contraction operator --------
class ContractGO(GraphOperator.GraphOperator):

    def __init__(self, n_vertices, n_loops, even_edges=True):

        directory = data_directory_even if even_edges else data_directory_odd
        s = "contract%d_%d.txt" % (n_vertices, n_loops)
        file_name = os.path.join(directory, s)

        domain = OrdinaryGVS(n_vertices, n_loops, even_edges=even_edges)
        target = OrdinaryGVS(n_vertices - 1, n_loops, even_edges=even_edges)

        super(ContractGO, self).__init__(file_name, domain, target)

    def _operate_on(self,G):
        image=[]
        for (i, e) in enumerate(G.edges(labels=False)):
            (u, v) = e
            # permute u,v to position 0,1
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

            pp = [j-1 for j in (Permutation([j+1 for j in p]).inverse())]
            sgn = self.domain.perm_sign(G, pp)
            G1 = copy(G)
            G1.relabel(pp, inplace=True)

            for (j, ee) in enumerate(G1.edges(labels=False)):
                a, b = ee
                G1.set_edge_label(a,b,j)
            # now delete the first edge and join the vertices 0 and 1
            G1.merge_vertices([0,1])
            G1.relabel(list(range(0,G1.order())), inplace = True)

            if not self.domain.even_edges:
                p = [j for (a, b, j) in G1.edges()]
                #print(p)
                #sgn *= Permutation(p).signature()

            image.append((G1, sgn))
        return image

    def get_work_estimate(self):
        # give estimate of number of graphs
        return self.domain.get_work_estimate() * self.domain.n_edges
