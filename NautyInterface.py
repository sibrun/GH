"""Provides Nauty interfaces to generate graphs.
For details, see: http://users.cecs.anu.edu.au/~bdm/nauty/
"""

__all__ = ['list_simple_graphs', 'list_bipartite_graphs']

from sage.all import *
import os
import tempfile
import Log

logger = Log.logger.getChild('nauty_interface')

def list_simple_graphs(n_vertices, n_edges, onlyonevi=True):
    """Creates a list of simple 1vi graphs with at least trivalent vertices.

    :param n_vertices: non-negative int: Number of vertices.
    :param n_edges: non-negative int: Number of edges.
    :param onlyonevi: TODO
    :return: list(sage.Graph): List of generated sage graphs.
    """
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    nauty_string = ("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (n_vertices, n_edges, n_edges)
    #logger.warn('call nauty to generate simple graphs: ' + nauty_string)
    graph_list = list(graphs.nauty_geng(nauty_string))
    #if len(graph_list) == 0:
        #print('Call nauty to generate bipartite graphs: ' + nauty_string)
        #print('List of bipartite graphs generated using nauty has zero length')
        #logger.warn('Call nauty to generate bipartite graphs: ' + nauty_string)
        #logger.warn('List of simple graphs generated using nauty has zero length')
    return graph_list


def list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges):
    """Creates a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
    vertices of the second colour have degree in the range min_deg_2:max_deg_2.
    :param n_vertices_1: non-negative int: Number of vertices of the first colour.
    :param n_vertices_2: non-negative int: Number of vertices of the second colour.
    :param deg_range_1: tuple(non-negative int, non-negative int): (min, max) degree of the vertices of the first colour.
    :param deg_range_2: tuple(non-negative int, non-negative int): (min, max) degree of the vertices of the second colour.
    :param n_edges: non-negative int: Number of edges of the bipartite graph.
    :return: list(sage.Graph): List of generated sage graphs.
    """
    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2
    # z switch prevents multiple hairs and multiple edges
    with tempfile.NamedTemporaryFile() as f:
        nauty_command = 'genbgL -czlq -d%d:%d -D%d:%d %d %d %d:%d %s' % \
                   (min_deg_1, min_deg_2, max_deg_1, max_deg_2, n_vertices_1, n_vertices_2, n_edges, n_edges, f.name)
        #logger.warn('call nauty to generate bipartite graphs: ' + nauty_command)
        os.system(nauty_command)
        list_g6 = f.read().splitlines()
        #if len(list_g6) == 0:
            #print('Call nauty to generate bipartite graphs: ' + nauty_command)
            #print('List of bipartite graphs generated using nauty has zero length')
            #logger.warn('Call nauty to generate bipartite graphs: ' + nauty_command)
            #logger.warn('List of bipartite graphs generated using nauty has zero length')
    return [Graph(g6) for g6 in list_g6]
