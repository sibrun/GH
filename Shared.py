"""This module contains shared code."""

__all__ = ['Perm', 'OrderedDict', 'enumerate_edges', 'edge_perm_sign', 'shifted_edge_perm_sign', 'permute_to_left',
           'matrix_norm', 'power_2']

from sage.all import *
import collections


class Perm:
    def __init__(self, p):
        self.p = p

    def inverse(self):
        inverse = [0] * len(self.p)
        for i, p in enumerate(self.p):
            inverse[p] = i
        return inverse
        #return [j-1 for j in Permutation([j+1 for j in self.p]).inverse()]

    def signature(self):
        return Permutation([j+1 for j in self.p]).signature()

    @classmethod
    def shifted(cls, p):
        pmin = min(p)
        return cls([j - pmin for j in p])


class OrderedDict(collections.OrderedDict):
    def __init__(self, *args):
        super(OrderedDict, self).__init__(*args)

    def __str__(self):
        s = '('
        for (key, value) in self.items():
            s += '%s: %s ' % (key, str(value))
        s += ')'
        return s

    def get_value_tuple(self):
        return tuple(self.values())


def enumerate_edges(graph):
    for (j, e) in enumerate(graph.edges(labels=False)):
        a, b = e
        graph.set_edge_label(a, b, j)


def edge_perm_sign(graph):
    p = [j for (a, b, j) in graph.edges()]
    return Perm(p).signature()


def shifted_edge_perm_sign(graph):
    p = [j for (a, b, j) in graph.edges()]
    return Perm.shifted(p).signature()


def permute_to_left((u, v), vertex_range):
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
    return sqrt(sum(map(power_2, M.list())))


def power_2(x):
    return x*x