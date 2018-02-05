from sage.all import *
import os
import StoreLoad as SL

temp_folder = 'temp'
temp_g6_file = 'temp.txt'


def get_temp_file_path():
    path = os.path.join(temp_folder, temp_g6_file)
    SL.generate_path(path)
    return path


'''creates a list of simple 1vi graphs with at least trivalent vertices'''
def list_simple_graphs(n_vertices, n_edges, onlyonevi=True):
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    nauty_string = ("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (n_vertices, n_edges, n_edges)
    return list(graphs.nauty_geng(nauty_string))


'''creates a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
vertices of the second colour have degree in the range min_deg_2:max_deg_2'''
def list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges):
    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2
    temp_file = get_temp_file_path()

    # z switch prevents multiple hairs and multiple edges
    nauty_command = 'genbg -czl -d%d:%d -D%d:%d %d %d %d:%d %s' % \
                   (min_deg_1, min_deg_2, max_deg_1, max_deg_2, n_vertices_1, n_vertices_2, n_edges, n_edges, temp_file)
    os.system(nauty_command)
    list_g6 = SL.load_string_list(temp_file)
    return [Graph(g6) for g6 in list_g6]



'''

function listG(nVertices, nLoops, onlyonevi=true)
  nEdges = nLoops + nVertices -1
  if (3*nVertices > 2*nEdges) || (nEdges > nVertices*(nVertices-1)/2 )
    # impossible
    return SmallGraph[]
  end
  tempFile = get_temp_file_name()
  run(`$geng $(onlyonevi?"-Cd3":"-cd3") $nVertices $nEdges:$nEdges $tempFile`)
  lll=readAllLines(tempFile)
  return [parse_graph6(l) for l in lll]
end


function parse_graphT(s)
    a = split(s)
    nEdges = parse(Int,a[2])
    nVert = parse(Int,a[1])
    G = small_multi_graph(nVert)
    for i = 1:nEdges
        u = parse(Int,a[3*i]) +1
        v = parse(Int,a[3*i+1]) +1
        emultiplicity = parse(Int,a[3 * i + 2])
        add_edge!(G,u,v, emultiplicity)
    end
    return G
end

function get_generating_graphs(self::HairyGraphVectorSpace)
    # Idea: produce all bipartite graphs, the second color being either of degree 1 or 2
    # degree 1 pieces are hairs, degree 2 vertices are edges and are removed later
    # z switch prevents multiple hairs and multiple edges

    nEdges = self.nLoops + self.nVertices -1
    nEdgesBiP = self.nHairs+2*nEdges
    tempFile = get_temp_file_name()
    run(`$genbg -czl -d3:1 -D$(nEdges+1):2 $(self.nVertices) $(self.nHairs+nEdges) $nEdgesBiP:$nEdgesBiP $tempFile`)
    lll=readAllLines(tempFile)
    return [bip_to_ordinary(self, parse_graph6(l)) for l in lll]
end
'''