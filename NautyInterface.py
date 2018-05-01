from sage.all import *
import os
import tempfile
import Log

logger = Log.logger.getChild('nauty_interface')

def list_simple_graphs(n_vertices, n_edges, onlyonevi=True):
    """creates a list of simple 1vi graphs with at least trivalent vertices"""
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    nauty_string = ("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (n_vertices, n_edges, n_edges)
    #logger.warn('call nauty to generate simple graphs: ' + nauty_string)
    graph_list = list(graphs.nauty_geng(nauty_string))
    if len(graph_list) == 0:
        logger.warn('list of simple graphs generated using nauty has zero length')
    return graph_list


def list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges):
    """creates a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
    vertices of the second colour have degree in the range min_deg_2:max_deg_2"""

    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2
    # z switch prevents multiple hairs and multiple edges
    with tempfile.NamedTemporaryFile() as f:
        nauty_command = 'genbg -czlq -d%d:%d -D%d:%d %d %d %d:%d %s' % \
                   (min_deg_1, min_deg_2, max_deg_1, max_deg_2, n_vertices_1, n_vertices_2, n_edges, n_edges, f.name)
        #logger.warn('call nauty to generate bipartite graphs: ' + nauty_command)
        os.system(nauty_command)
        list_g6 = f.read().splitlines()
        if len(list_g6) == 0:
            logger.warn('list of bipartite graphs generated using nauty has zero length')
    return [Graph(g6) for g6 in list_g6]
