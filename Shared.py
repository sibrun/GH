import os
from sage.all import *


class Perm():
    def __init__(self, p):
        self.p = p

    def inverse(self):
        return [j-1 for j in Permutation([j+1 for j in self.p]).inverse()]

    def sign(self):
        return Permutation([j+1 for j in self.p]).signature()


class NotBuiltError(RuntimeError):
    pass

class RefError(RuntimeError):
    pass

def get_log_path(log_dir, log_file):
    log_path = os.path.join(log_dir, log_file)
    log_dir = os.path.dirname(log_path)
    if not os.path.exists(log_dir):
        os.makedirs(log_dir)
    return log_path