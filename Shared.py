from sage.all import *


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

