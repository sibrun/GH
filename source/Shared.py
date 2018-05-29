"""Provide shared code."""

__all__ = ['Perm', 'OrderedDict', 'enumerate_edges', 'edge_perm_sign', 'shifted_edge_perm_sign', 'permute_to_left',
           'matrix_norm', 'power_2']

from sage.all import *
import collections


class Perm:
    """Permutation starting at zero."""

    def __init__(self, p):
        """Initialize the permutation.

        :param p: list(non-negative int): Image of the permutation with consecutive indices starting at zero.
        """
        self.p = p

    def inverse(self):
        """Return the inverse permutation.

        :return: list(non-negative int): Image of the inverse permutation.
        """
        inverse = [0] * len(self.p)
        for i, p in enumerate(self.p):
            inverse[p] = i
        return inverse
        #return [j-1 for j in Permutation([j+1 for j in self.p]).inverse()]

    def signature(self):
        """Returns the sign of the permutation.

        :return: int: Sign of the permutation.
        """
        return Permutation([j+1 for j in self.p]).signature()

    @classmethod
    def shifted(cls, p):
        """Generates a permutation not necessarily starting at zero.

        Caution: The inverse of the returned permutation does not correspond to the inverse of the permutation induced by
        p but rather to the corresponding inverse permutation of the shifted indices starting at zero.

        :param p: list(non-negative int): Image of the permutation with consecutive indices.
        :return: Perm: Permutation of consecutive indices starting at zero.
        """
        pmin = min(p)
        return cls([j - pmin for j in p])


class OrderedDict(collections.OrderedDict):
    """Ordered dictionary.

    """
    def __init__(self, *args):
        super(OrderedDict, self).__init__(*args)

    def __str__(self):
        """Returns the string of the ordered dictionary.

        :return: str: Ordered dictionary as string.
        """
        s = '('
        for (key, value) in self.items():
            s += '%s: %s ' % (key, str(value))
        s += ')'
        return s

    def get_value_tuple(self):
        """Returns the values as tuple.

        :return: tuple: Values of the ordered dictionary as tuple.
        """
        return tuple(self.values())


def enumerate_edges(graph):
    """Labels the edges of the graph lexicographically.

    :param graph: sage.Graph: Input graph.
    """
    for (j, e) in enumerate(graph.edges(labels=False)):
        a, b = e
        graph.set_edge_label(a, b, j)


def edge_perm_sign(graph):
    """Returns the sign of the permutation induced by the order of edges of the graph.

    Edge labels are supposed to start at zero.

    :param graph: sage.Graph: Input graph.
    :return: int: Sign of the permutation induced by the order of edges of the graph.
    """
    p = [j for (a, b, j) in graph.edges()]
    return Perm(p).signature()


def shifted_edge_perm_sign(graph):
    """Returns the sign of the permutation induced by the order of edges of the graph.

    Edge labels don't have to start at zero.

    :param graph: sage.Graph: Input graph.
    :return: int: Sign of the permutation induced by the order of edges of the graph.
    """
    p = [j for (a, b, j) in graph.edges()]
    return Perm.shifted(p).signature()


def permute_to_left(pair, vertex_range):
    """Permutes u,v to the left of the range vertex_range and returns the induced permutation.

    :param vertex_range: range: Range of vertices.
    :return: list(non-negative int): Permutation which permutes u,v to the left of the range vertex_range.
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
    """Returns the Frobenius norm of the matrix M.

    :param M: sage.Matrix: Input matrix.
    :return: non-negative int: Frobenius norm of the matrix M.
    """
    return sqrt(sum(map(power_2, M.list())))


def power_2(x):
    return x*x