import unittest
import OrdinaryGraphComplex as OGC
from sage.all import *

reload(OGC)


class OGCTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_perm_sign(self):
        ogc_even = OGC.OrdinaryGraphVectorSpace(6, 5, True)
        ogc_odd = OGC.OrdinaryGraphVectorSpace(6, 5, False)
        G = graphs.WheelGraph(5)
        p = [0, 2, 3, 4, 5, 1]
        self.assertTrue(ogc_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        self.assertTrue(ogc_even.perm_sign(G, p) == 1, 'incorrect permutation sign')
        p = [0, 1, 5, 4, 3, 2]
        self.assertTrue(ogc_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        self.assertTrue(ogc_even.perm_sign(G, p) == -1, 'incorrect permutation sign')
        p = [0, 1, 2, 4, 3, 5]
        self.assertTrue(ogc_odd.perm_sign(G, p) == -1, 'incorrect permutation sign')
        self.assertTrue(ogc_even.perm_sign(G, p) == 1, 'incorrect permutation sign')


def suite():
    suite = unittest.TestSuite()
    suite.addTest(OGCTestCase('test_perm_sign'))
    return suite


if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())