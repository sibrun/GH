"""Provide Nauty interfaces to generate graphs.
For details, see: http://users.cecs.anu.edu.au/~bdm/nauty/
"""

__all__ = ['list_simple_graphs', 'list_bipartite_graphs']

from sage.all import *
import os
import tempfile
import Log
import Parameters
import StoreLoad

logger = Log.logger.getChild('nauty_interface')


def run_sys_cmd(cmd):
    """Run system command, raising an error if an error exit code is detected."""
    ret = os.system(cmd)
    if ret != 0:
        raise RuntimeError(f"Nonzero exit code from: {cmd}")


def list_simple_graphs(n_vertices, n_edges, onlyonevi=True):
    """Create a list of simple 1vi graphs with at least trivalent vertices.

    :param n_vertices: Number of vertices.
    :type n_vertices: int
    :param n_edges: Number of edges.
    :type n_edges: int
    :param onlyonevi: TODO
    :return: List of generated sage graphs.
    :rtype: list(Graph)
    """
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    nauty_string = ("-Cd3" if onlyonevi else "-cd3") + \
        " %d %d:%d" % (n_vertices, n_edges, n_edges)
    #logger.warning('call nauty to generate simple graphs: ' + nauty_string)
    graph_list = graphs.nauty_geng(nauty_string)
    # graph_list = list(graphs.nauty_geng(nauty_string))
    # if len(graph_list) == 0:
    #print('Call nauty to generate bipartite graphs: ' + nauty_string)
    #print('List of bipartite graphs generated using nauty has zero length')
    #logger.warning('Call nauty to generate bipartite graphs: ' + nauty_string)
    #logger.warning('List of simple graphs generated using nauty has zero length')
    return graph_list

def list_simple_graphs_1(n_vertices, n_edges, onlyonevi=True):
    """Create a list of simple 1vi graphs with at least 1-valent vertices.

    :param n_vertices: Number of vertices.
    :type n_vertices: int
    :param n_edges: Number of edges.
    :type n_edges: int
    :param onlyonevi: TODO
    :return: List of generated sage graphs.
    :rtype: list(Graph)
    """
    if n_vertices <= 0 or n_edges <= 0 or n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    nauty_string = ("-Cd1" if onlyonevi else "-cd1") + \
        " %d %d:%d" % (n_vertices, n_edges, n_edges)
    graph_list = graphs.nauty_geng(nauty_string)
    return graph_list if graph_list else []


def list_simple_graphs_valence(n_vertices, n_edges, max_valence, onlyonevi=True):
    """Create a list of simple 1vi graphs with at least trivalent vertices.

    :param n_vertices: Number of vertices.
    :type n_vertices: int
    :param n_edges: Number of edges.
    :type n_edges: int
    :param onlyonevi: TODO
    :return: List of generated sage graphs.
    :rtype: list(Graph)
    """
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    nauty_string = ("-Cd3" if onlyonevi else "-cd3") + (f" -D{max_valence}") + \
        " %d %d:%d" % (n_vertices, n_edges, n_edges)
    #logger.warning('call nauty to generate simple graphs: ' + nauty_string)
    graph_list = graphs.nauty_geng(nauty_string)
    # graph_list = list(graphs.nauty_geng(nauty_string))
    # if len(graph_list) == 0:
    #print('Call nauty to generate bipartite graphs: ' + nauty_string)
    #print('List of bipartite graphs generated using nauty has zero length')
    #logger.warning('Call nauty to generate bipartite graphs: ' + nauty_string)
    #logger.warning('List of simple graphs generated using nauty has zero length')
    return graph_list


def list_bipartite_graphs(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges):
    """Create a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
    vertices of the second colour have degree in the range min_deg_2:max_deg_2.

    :param n_vertices_1: Number of vertices of the first colour.
    :type n_vertices_1: int
    :param n_vertices_2: Number of vertices of the second colour.
    :type n_vertices_2: int
    :param deg_range_1: (min, max) degree of the vertices of the first colour.
    :type deg_range_1: tuple(int, int)
    :param deg_range_2: (min, max) degree of the vertices of the second colour.
    :type deg_range_2: tuple(int, int)
    :param n_edges: Number of edges of the bipartite graph.
    :type n_edges: int
    :return: List of generated sage graphs.
    :rtype: list(Graph)
    """
    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2
    # z switch prevents multiple hairs and multiple edges
    StoreLoad.makedirs(Parameters.temp_folder)
    with tempfile.NamedTemporaryFile(dir=Parameters.temp_folder) as f:
        nauty_command = 'genbgL -czlq -d%d:%d -D%d:%d %d %d %d:%d %s' % \
            (min_deg_1, min_deg_2, max_deg_1, max_deg_2,
             n_vertices_1, n_vertices_2, n_edges, n_edges, f.name)
        #logger.warning('call nauty to generate bipartite graphs: ' + nauty_command)
        run_sys_cmd(nauty_command)
        txt = f.read()
        if type(txt) is not str:
            txt = txt.decode("ascii")
        list_g6 = txt.splitlines()
        print(
            f"Call to nauty: {nauty_command} generated {len(list_g6)} graphs.")
        #list_g6 = f.read().splitlines()
        # if len(list_g6) == 0:
        #print('Call nauty to generate bipartite graphs: ' + nauty_command)
        #print('List of bipartite graphs generated using nauty has zero length')
        #logger.warning('Call nauty to generate bipartite graphs: ' + nauty_command)
        #logger.warning('List of bipartite graphs generated using nauty has zero length')
    return (Graph(g6) for g6 in list_g6)


def list_bipartite_graphs2(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges):
    """
    Same as before, but with at most one common neighbor ...
    Create a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
    vertices of the second colour have degree in the range min_deg_2:max_deg_2.

    :param n_vertices_1: Number of vertices of the first colour.
    :type n_vertices_1: int
    :param n_vertices_2: Number of vertices of the second colour.
    :type n_vertices_2: int
    :param deg_range_1: (min, max) degree of the vertices of the first colour.
    :type deg_range_1: tuple(int, int)
    :param deg_range_2: (min, max) degree of the vertices of the second colour.
    :type deg_range_2: tuple(int, int)
    :param n_edges: Number of edges of the bipartite graph.
    :type n_edges: int
    :param maxneighbors: Maximum number of common neighbors of second color vertices
    :return: List of generated sage graphs.
    :rtype: list(Graph)
    """
    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2
    # z switch prevents multiple hairs and multiple edges
    StoreLoad.makedirs(Parameters.temp_folder)
    with tempfile.NamedTemporaryFile(dir=Parameters.temp_folder) as f:
        nauty_command = 'genbgL -clq -Z1 -d%d:%d -D%d:%d %d %d %d:%d %s' % \
            (min_deg_1, min_deg_2, max_deg_1, max_deg_2,
             n_vertices_1, n_vertices_2, n_edges, n_edges, f.name)
        # print(nauty_command)
        #logger.warning('call nauty to generate bipartite graphs: ' + nauty_command)
        run_sys_cmd(nauty_command)
        txt = f.read()
        if type(txt) is not str:
            txt = txt.decode("ascii")
        list_g6 = txt.splitlines()
        print(
            f"Call to nauty: {nauty_command} generated {len(list_g6)} graphs.")
        # list_g6 = f.read().decode("utf-8").splitlines()
        # print(list_g6)
        # if len(list_g6) == 0:
        #print('Call nauty to generate bipartite graphs: ' + nauty_command)
        #print('List of bipartite graphs generated using nauty has zero length')
        #logger.warning('Call nauty to generate bipartite graphs: ' + nauty_command)
        #logger.warning('List of bipartite graphs generated using nauty has zero length')
    return (Graph(g6) for g6 in list_g6)


def list_bipartite_graphs3(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges, n_maxneighbors):
    """
    Same as before, but with  n_maxcommon_neighbors many common neighbors ...
    Create a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
    vertices of the second colour have degree in the range min_deg_2:max_deg_2.

    :param n_vertices_1: Number of vertices of the first colour.
    :type n_vertices_1: int
    :param n_vertices_2: Number of vertices of the second colour.
    :type n_vertices_2: int
    :param deg_range_1: (min, max) degree of the vertices of the first colour.
    :type deg_range_1: tuple(int, int)
    :param deg_range_2: (min, max) degree of the vertices of the second colour.
    :type deg_range_2: tuple(int, int)
    :param n_edges: Number of edges of the bipartite graph.
    :type n_edges: int
    :param n_maxneighbors: Maximum number of common neighbors of second color vertices
    :return: List of generated sage graphs.
    :rtype: list(Graph)
    """
    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2
    # z switch prevents multiple hairs and multiple edges

    StoreLoad.makedirs(Parameters.temp_folder)
    with tempfile.NamedTemporaryFile(dir=Parameters.temp_folder) as f:
        nauty_command = 'genbgL -clq -Z%d -d%d:%d -D%d:%d %d %d %d:%d %s' % \
            (n_maxneighbors, min_deg_1, min_deg_2, max_deg_1, max_deg_2,
             n_vertices_1, n_vertices_2, n_edges, n_edges, f.name)
        # print(nauty_command)
        #logger.warning('call nauty to generate bipartite graphs: ' + nauty_command)
        run_sys_cmd(nauty_command)
        txt = f.read()
        if type(txt) is not str:
            txt = txt.decode("ascii")
        list_g6 = txt.splitlines()
        print(
            f"Call to nauty: {nauty_command} generated {len(list_g6)} graphs.")
        # list_g6 = f.read().decode("utf-8").splitlines()
        # print(list_g6)
        # if len(list_g6) == 0:
        #print('Call nauty to generate bipartite graphs: ' + nauty_command)
        #print('List of bipartite graphs generated using nauty has zero length')
        #logger.warning('Call nauty to generate bipartite graphs: ' + nauty_command)
        #logger.warning('List of bipartite graphs generated using nauty has zero length')
    return (Graph(g6) for g6 in list_g6)


def list_bipartite_graphs_disc(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges, n_maxneighbors):
    """
    Same as before, but including disconnected graphs ...
    Create a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
    vertices of the second colour have degree in the range min_deg_2:max_deg_2.

    :param n_vertices_1: Number of vertices of the first colour.
    :type n_vertices_1: int
    :param n_vertices_2: Number of vertices of the second colour.
    :type n_vertices_2: int
    :param deg_range_1: (min, max) degree of the vertices of the first colour.
    :type deg_range_1: tuple(int, int)
    :param deg_range_2: (min, max) degree of the vertices of the second colour.
    :type deg_range_2: tuple(int, int)
    :param n_edges: Number of edges of the bipartite graph.
    :type n_edges: int
    :param n_maxneighbors: Maximum number of common neighbors of second color vertices
    :return: List of generated sage graphs.
    :rtype: list(Graph)
    """
    print("list_bipartite_graphs_disc", n_vertices_1, n_vertices_2,
          deg_range_1, deg_range_2, n_edges, n_maxneighbors)
    (min_deg_1, max_deg_1) = deg_range_1
    (min_deg_2, max_deg_2) = deg_range_2
    # z switch prevents multiple hairs and multiple edges

    StoreLoad.makedirs(Parameters.temp_folder)
    with tempfile.NamedTemporaryFile(dir=Parameters.temp_folder) as f:
        nauty_command = 'genbgL -lq -Z%d -d%d:%d -D%d:%d %d %d %d:%d %s' % \
            (n_maxneighbors, min_deg_1, min_deg_2, max_deg_1, max_deg_2,
             n_vertices_1, n_vertices_2, n_edges, n_edges, f.name)
        # print(nauty_command)
        #logger.warning('call nauty to generate bipartite graphs: ' + nauty_command)
        run_sys_cmd(nauty_command)
        txt = f.read()
        if type(txt) is not str:
            txt = txt.decode("ascii")
        list_g6 = txt.splitlines()
        print(
            f"Call to nauty: {nauty_command} generated {len(list_g6)} graphs.")
        # list_g6 = f.read().decode("utf-8").splitlines()
        # print(list_g6)
        # if len(list_g6) == 0:
        #print('Call nauty to generate bipartite graphs: ' + nauty_command)
        #print('List of bipartite graphs generated using nauty has zero length')
        #logger.warning('Call nauty to generate bipartite graphs: ' + nauty_command)
        #logger.warning('List of bipartite graphs generated using nauty has zero length')
    return (Graph(g6) for g6 in list_g6)


# def list_bipartite_graphs3(n_vertices_1, n_vertices_2, deg_range_1, deg_range_2, n_edges):
#     """
#     Same as before, but with at most one common neighbor ...
#     Create a list of bipartite graphs, vertices of the first colour have degree in the range min_deg_1:max_deg_1,
#     vertices of the second colour have degree in the range min_deg_2:max_deg_2.

#     :param n_vertices_1: Number of vertices of the first colour.
#     :type n_vertices_1: int
#     :param n_vertices_2: Number of vertices of the second colour.
#     :type n_vertices_2: int
#     :param deg_range_1: (min, max) degree of the vertices of the first colour.
#     :type deg_range_1: tuple(int, int)
#     :param deg_range_2: (min, max) degree of the vertices of the second colour.
#     :type deg_range_2: tuple(int, int)
#     :param n_edges: Number of edges of the bipartite graph.
#     :type n_edges: int
#     :param maxneighbors: Maximum number of common neighbors of second color vertices
#     :return: List of generated sage graphs.
#     :rtype: list(Graph)
#     """
#     (min_deg_1, max_deg_1) = deg_range_1
#     (min_deg_2, max_deg_2) = deg_range_2
#     # z switch prevents multiple hairs and multiple edges
#     with tempfile.NamedTemporaryFile() as f:
#         nauty_command = 'genbgL -lq -Z1 -d%d:%d -D%d:%d %d %d %d:%d %s' % \
#                    (min_deg_1, min_deg_2, max_deg_1, max_deg_2, n_vertices_1, n_vertices_2, n_edges, n_edges, f.name)
#         #print(nauty_command)
#         #logger.warning('call nauty to generate bipartite graphs: ' + nauty_command)
#         os.system(nauty_command)
#         list_g6 = f.read().decode("utf-8").splitlines()
#         #print(list_g6)
#         #if len(list_g6) == 0:
#             #print('Call nauty to generate bipartite graphs: ' + nauty_command)
#             #print('List of bipartite graphs generated using nauty has zero length')
#             #logger.warning('Call nauty to generate bipartite graphs: ' + nauty_command)
#             #logger.warning('List of bipartite graphs generated using nauty has zero length')
#     return [Graph(g6) for g6 in list_g6]
