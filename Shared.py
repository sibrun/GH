import os
import pickle
import scipy.sparse as sparse
from sage.all import *


class Perm:
    def __init__(self, p):
        self.p = p

    def inverse(self):
        return [j-1 for j in Permutation([j+1 for j in self.p]).inverse()]

    def sign(self):
        return Permutation([j+1 for j in self.p]).signature()

def nauty_geng(self, n_vertices, n_edges, onlyonevi=True):
    if n_vertices <= 0 or n_edges <= 0 or 3 * n_vertices > 2 * n_edges or  self.n_edges > self.n_vertices * (self.n_vertices - 1) / 2:
        return []
    return list(graphs.nauty_geng(("-Cd3" if onlyonevi else "-cd3") + " %d %d:%d" % (n_vertices, n_edges, n_edges)))


class NotBuiltError(RuntimeError):
    pass


class RefError(RuntimeError):
    pass


class FileNotExistingError(RuntimeError):
    pass


def get_path_from_current(*paths):
    return os.path.join(os.getcwd(), *paths)


def generate_path(path):
    directory = os.path.dirname(path)
    if not os.path.exists(directory):
        os.makedirs(directory)


def store_string_list(L, path, header=None):
    generate_path(path)
    with open(path,'w') as f:
        if header is not None:
            f.write(header + '\n')
        for x in L:
            f.write(x + '\n')


def load_string_list(path, header=True):
    if not os.path.exists(path):
        raise FileNotExistingError("Cannot load from %s: The file does not exist" % str(path))
    with open(path, 'r') as f:
        if header:
            H = f.readline()
        L = f.read().splitlines()
    if header:
        return (H, L)
    else:
        return L


def load_header(path):
    if not os.path.exists(path):
        raise FileNotExistingError("Cannot load from %s: The file does not exist" % str(path))
    with open(path, 'r') as f:
        return f.readline()


def store_list_of_header_lists(LHL, path):
    generate_path(path)
    with open(path, 'w') as f:
        for HL in LHL:
            (H,L) = HL
            if H is not None:
                f.write(H + '\n')
            for x in L:
                f.write(x + '\n')


def pickle_store(Ob, path):
    generate_path(path)
    with open(path,'wb') as f:
        pickle.dump(Ob, f)


def pickle_load(path):
    if not os.path.exists(path):
        raise FileNotExistingError("Cannot load from %s: The file does not exist" % str(path))
    with open(path, 'rb') as f:
        Ob = pickle.load(f)
    return Ob


def sparse_inverse(M):
    (m ,n) = M.shape
    if m != n:
        raise ValueError("Cannot compute the inverse: Not a square matrix")
    return sparse.csc_matrix(sparse.linalg.spsolve(M, sparse.eye(m))).astype(int)