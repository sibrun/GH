# Temporarily
import os


ordinaryDataDirOdd = DATA_DIR + "/ordinarydata/oddedge/"
ordinaryDataDirEven = DATA_DIR + "/ordinarydata/evenedge/"
imgBaseDir = "img/"

#- ---auxiliary functions (to be moved elsewhere)------------
def listG(nVertices, nLoops, onlyonevi=True):
    """
    creates a list of simple 1vi graphs with at least trivalent vertices
    """
    nEdges = nLoops + nVertices -1
    if (3*nVertices > 2*nEdges) || (nEdges > nVertices*(nVertices-1)/2 ):
        # impossible
        return []

    gL = list(graphs.nauty_geng(("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (nVertices, nEdges,nEdges)))
    return gL

#-----------------------

class OrdinaryGraphVectorSpace(GraphVectorSpace):
    def __init__(self, nVertices, nLoops, evenEdges):
        self.nVertices = nVertices
        self.nLoops = nLoops
        self.evenEdges = evenEdges


    def get_file_name(self):
        dataDir = ordinaryDataDirEven if self.evenEdges else ordinaryDataDirOdd
        s = "gra" + self.nVertices + "_" + self.nLoops + ".g6"
        return os.path.join(dataDir, s)


    def get_svg_dir(self):
        dataDir = ordinaryDataDirEven if self.evenEdges else ordinaryDataDirOdd
        s = "imgs%d_%d/" % (self.nVertices, self.nLoops)
        return os.path.join(dataDir, imgBaseDir, s)

    def get_color_counts(self):
        return None # no coloring


    def is_valid(self):
        nEdges = self.nLoops + self.nVertices -1
        # at least trivalent, and simple
        return  (3*self.nVertices <= 2*nEdges) and self.nVertices > 0 and self.nLoops >= 0 and nEdges <= self.nVertices*(self.nVertices-1)/2


    def get_generating_graphs(self):
        # generate List of unmarked graphs
        #println(self.nVertices)
        LL = listG(self.nVertices, self.nLoops)
        #LL = listG(self.nVertices, self.nLoops, false)
        #LL=[]
        return LL



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
            pp = range(nEdges)
            G2 = permuteGraph(G,p)
            edge_idx2 = -ones(Int64, nVert, nVert)
            for (j, e) = enumerate(edges(G2))
                u,v = e
                edge_idx2[u,v] = j
                edge_idx2[v,u] = j
            end
            for (j, e) = enumerate(edges(G))
                u,v = e
                pp[j] = edge_idx2[p[u],p[v]]
            end

            #println(pp)
            return permSign(pp)

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

