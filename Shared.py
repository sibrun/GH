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