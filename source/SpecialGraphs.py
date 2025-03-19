"""Contains routines that generate some common types of graphs.
"""

import OrdinaryGraphComplex
import CHairyGraphComplex
from sage.all import Graph


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
    p is the zero-based permutation connecting both halves of the graph.
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


def forested_morita_tetrahedron(ps):
    """Generates the graphs whose linear combinations conjecturally spans
    the degree 8 cohomology of the 7-loop odd_e forested graph complex.
    They have 12 vertices and 10 unmarked edges and 8 marked edges.
    ps are 4 zero-based permutations of (0,1,2) stating how the edges should be connected.
    Mind that unmarked edges are represented by bivalent vertices in the forested
    graph complex.
    """
    n_unmarked = 10
    n_vertices = 12
    # n_marked = 8
    G = Graph(n_vertices + n_unmarked)

    for j in range(4):
        for jj in range(2):
            # add marked edges
            G.add_edge(3*j+jj, 3*j+jj+1)

    # unmarked edges in the four Lie-vertices of the tetrahedron
    G.add_edge(0, n_vertices)
    G.add_edge(2, n_vertices)
    G.add_edge(3, n_vertices+1)
    G.add_edge(5, n_vertices+1)
    G.add_edge(6, n_vertices+2)
    G.add_edge(8, n_vertices+2)
    G.add_edge(9, n_vertices+3)
    G.add_edge(11, n_vertices+3)

    # unmarked edges between Lie vertices
    unmarked_ind = n_vertices + 4
    lie_inds = [0, 0, 0, 0]
    for i in range(4):
        for j in range(i+1, 4):
            G.add_edge(ps[i][lie_inds[i]] + i*3, unmarked_ind)
            G.add_edge(ps[j][lie_inds[j]] + j*3, unmarked_ind)
            unmarked_ind += 1
            lie_inds[i] = lie_inds[i]+1
            lie_inds[j] = lie_inds[j]+1

    return G
