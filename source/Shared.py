"""Provide shared code."""

__all__ = ['Perm', 'OrderedDict', 'enumerate_edges', 'edge_perm_sign', 'shifted_edge_perm_sign', 'permute_to_left',
           'matrix_norm', 'power_2']

from sage.all import *
import collections


class Perm:
    """Permutation starting at zero."""

    def __init__(self, p):
        """Initialize the permutation.

        :param p: Image of the permutation with consecutive indices starting at zero.
        :type p: list(int)
        """
        self.p = p

    def inverse(self):
        """Return the inverse permutation.

        :return: Image of the inverse permutation.
        :rtype: list(int)
        """
        inverse = [0] * len(self.p)
        for i, p in enumerate(self.p):
            inverse[p] = i
        return inverse
        #return [j-1 for j in Permutation([j+1 for j in self.p]).inverse()]

    def signature(self):
        """Return the sign of the permutation.

        :return: Sign of the permutation.
        :rtype: int
        """
        return Permutation([j+1 for j in self.p]).signature()

    @classmethod
    def shifted(cls, p):
        """Generate a permutation not necessarily starting at zero.

        :Note: The inverse of the returned permutation does not correspond to the inverse of the permutation induced by
               p but rather to the corresponding inverse permutation of the shifted indices starting at zero.

        :param p: Image of the permutation with consecutive indices.
        :type p: list(int)
        :return: Permutation of consecutive indices starting at zero.
        :rtype: Perm
        """
        pmin = min(p)
        return cls([j - pmin for j in p])


class OrderedDict(collections.OrderedDict):
    """Ordered dictionary.

    """
    def __init__(self, *args):
        super().__init__(*args)

    def __str__(self):
        """Return the string of the ordered dictionary.

        :return: Ordered dictionary as string.
        :rtype: str
        """
        s = '('
        for (key, value) in self.items():
            s += '%s: %s ' % (key, str(value))
        s += ')'
        return s

    def get_value_tuple(self):
        """Return the values as tuple.

        :return: Values of the ordered dictionary as tuple.
        :rtype: tuple
        """
        return tuple(self.values())


def enumerate_edges(graph):
    """Label the edges of the graph lexicographically.

    :param graph: Input graph.
    :type graph: Graph
    """
    for (j, e) in enumerate(graph.edges(labels=False)):
        a, b = e
        graph.set_edge_label(a, b, j)


def edge_perm_sign(graph):
    """Return the sign of the permutation induced by the order of edges of the graph.

    Edge labels are supposed to start at zero.

    :param graph: Input graph.
    :type graph: Graph
    :return: Sign of the permutation induced by the order of edges of the graph.
    :rtype: int
    """
    p = [j for (a, b, j) in graph.edges()]
    return Perm(p).signature()


def shifted_edge_perm_sign(graph):
    """Return the sign of the permutation induced by the order of edges of the graph.

    Edge labels don't have to start at zero.

    :param graph: Input graph.
    :type graph: Graph
    :return: Sign of the permutation induced by the order of edges of the graph.
    :rtype: int
    """
    p = [j for (a, b, j) in graph.edges()]
    # print("Edge perm", p)
    return Perm.shifted(p).signature()

def shifted_edge_perm_sign2(graph):
    """Return the sign of the permutation induced by the order of edges of the graph,
    that is the sign of the permutation needed to bring the edge labels in ascending order.
    :param graph: Input graph.
    :type graph: Graph
    :return: Sign of the permutation induced by the order of edges of the graph.
    :rtype: int
    """

    L = [(j,i) for i, (a, b, j) in enumerate(graph.edges())]
    L.sort()
    p = [i for j,i in L]
    return Perm.shifted(p).signature()

def permute_to_left(pair, vertex_range):
    """Permute pair to the left of the range vertex_range and returns the induced permutation.

    :param pair: Pair of vertex indices to be permuted to the left.
    :type pair: tuple(int, int)
    :param vertex_range: Range of vertices.
    :type vertex_range: range
    :return: Permutation which permutes u,v to the left of the range vertex_range.
    :rtype: list(int)
    """
    (u, v) = pair
    p = list(vertex_range)
    min_index = min(vertex_range)
    p[min_index] = u
    p[min_index + 1] = v
    idx = 2
    for j in vertex_range:
        if j == u or j == v:
            continue
        else:
            p[idx] = j
            idx += 1
    return Perm(p).inverse()


def matrix_norm(M):
    """Return the Frobenius norm of the matrix M.

    :param M: Input matrix.
    :type M: Matrix
    :return: Frobenius norm of the matrix M.
    :rtype: int
    """
    return sqrt(sum(map(power_2, M.list())))


def power_2(x):
    return x*x
