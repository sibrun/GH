"""Contains routines that generate some common types of graphs.
"""

import OrdinaryGraphComplex
import CHairyGraphComplex
from sage.all import *


def wheel_graph(nspokes):
    """Generates the wheel graphs, i.e., a central vertex connected to a ring of n vertices.
    The number of vertices will be nspokes+1. The central vertex is the first in the ordering.

        :param nspokes: The number of spokes in the wheel. 
        :type nspokes: int
    """
    G = Graph(nspokes+1)
    for j in range(nspokes):
        G.add_edge(0, j+1)
        G.add_edge(j+1, (j+2) % nspokes + 1)
    return G


def hedgehog_graph(nvertices):
    """Generates the hedgehog graphs, i.e., a ring of nvertices vertices, 
    each carrying one hair, with a unary vertex at the end.
    The trivalent vertices will be the first in the ordering

        :param nvertices: The number of trivalent vertices in the ring. 
        :type nvertices: int
    """
    G = Graph(2*nvertices)
    for j in range(nvertices):
        G.add_edge(j, j+nvertices)
        G.add_edge(j, (j+1) % nvertices)
    return G

def forested_ring_graph(n_marked_edges : int):
    """Generates the ring graphs that span the low degree cohomology 
    of the even_e forested graph complex.
    They have 2*n_marked_edges vertices and 2*n_marked_edges unmarked edges.
    Mind that unmarked edges are represented by bivalent vertices in the forested 
    graph complex.
    """
    n_unmarked = 2*n_marked_edges
    n_vertices = 2*n_marked_edges
    G = Graph(n_vertices + n_unmarked)

    for j in range(n_marked_edges):
        # add marked edges
        G.add_edge(2*j, 2*j+1)
        # unmarked edges
        G.add_edge(2*j, n_vertices + 2*j)
        G.add_edge(2*j+1, n_vertices + 2*j)
        G.add_edge(2*j+1, n_vertices + 2*j+1)
        G.add_edge((2*j+2) % n_vertices, n_vertices + 2*j+1)

    return G

def forested_morita_graph(k: int, p):
    """Generates the Morita graphs whose linear combinations span the low degree cohomology 
    of the odd_e forested graph complex.
    They have 2*k vertices and k+2 unmarked edges and 2*k-2 marked edges.
    p is the zero-based permutation connecting both halfs of the graph.
    Mind that unmarked edges are represented by bivalent vertices in the forested 
    graph complex.
    """
    n_unmarked = k+2
    n_vertices = 2*k
    # n_marked = 2*k-2
    G = Graph(n_vertices + n_unmarked)

    for j in range(k-1):
        # add marked edges
        G.add_edge(j, j+1)
        G.add_edge(k+j, k+j+1)

    # unmarked edges
    G.add_edge(0, n_vertices)
    G.add_edge(k-1, n_vertices)
    G.add_edge(k, n_vertices+1)
    G.add_edge(2*k-1, n_vertices+1)

    for j in range(k):
        G.add_edge(p[j], n_vertices + 2 + j)
        G.add_edge(k+j, n_vertices + 2 + j)

    return G