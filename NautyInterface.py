from sage.all import *
import os
import StoreLoad as SL

temp_folder = 'temp'


'''creates a list of simple 1vi graphs with at least trivalent vertices'''
def list_simple_graphs(n_vertices, n_edges, onlyonevi=True):
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    nauty_string = ("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (n_vertices, n_edges, n_edges)
    return list(graphs.nauty_geng(nauty_string))


'''creates a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
vertices of the second colour have degree in the range min_deg_2:max_deg_2'''
def list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges):
    s = 'bipartite_%d_%d_%d.txt' % (n_vertices_1, n_vertices_2, n_edges)
    temp_file = os.path.join(temp_folder, s)
    SL.generate_path(temp_file)

    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2

    # z switch prevents multiple hairs and multiple edges
    nauty_command = 'genbg -czlq -d%d:%d -D%d:%d %d %d %d:%d %s' % \
                   (min_deg_1, min_deg_2, max_deg_1, max_deg_2, n_vertices_1, n_vertices_2, n_edges, n_edges, temp_file)
    os.system(nauty_command)

    list_g6 = SL.load_string_list(temp_file)
    os.remove(temp_file)
    return [Graph(g6) for g6 in list_g6]
