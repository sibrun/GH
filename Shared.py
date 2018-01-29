from sage.all import *


class Perm:
    def __init__(self, p):
        self.p = p

    def inverse(self):
        return [j-1 for j in Permutation([j+1 for j in self.p]).inverse()]

    def sign(self):
        return Permutation([j+1 for j in self.p]).signature()


def nauty_geng(n_vertices, n_edges, onlyonevi=True):
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or n_edges > n_vertices * (n_vertices - 1) / 2:
        return []
    return list(graphs.nauty_geng(("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (n_vertices, n_edges, n_edges)))


class NotBuiltError(RuntimeError):
    pass


class RefError(RuntimeError):
    pass