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

    def sign(self):
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


def enumerate_edges(graph):
    for (j, e) in enumerate(graph.edges(labels=False)):
        a, b = e
        graph.set_edge_label(a, b, j)


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