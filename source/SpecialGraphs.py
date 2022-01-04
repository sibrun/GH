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
