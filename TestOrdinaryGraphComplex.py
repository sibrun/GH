import unittest
import OrdinaryGraphComplex as OGC
import TestGraphVectorSpace as TGVS
from sage.all import *

reload(OGC)


class OGCTestCase(unittest.TestCase):

    def test_perm_sign(self):

        ogc_even = OGC.OrdinaryGraphVectorSpace(6, 5, evenEdges=True)
        ogc_odd = OGC.OrdinaryGraphVectorSpace(6, 5, evenEdges=False)
        G = graphs.WheelGraph(5)
        p = [0, 2, 3, 4, 5, 1]
        self.assertTrue(ogc_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        self.assertTrue(ogc_even.perm_sign(G, p) == 1, 'incorrect permutation sign')
        p = [0, 1, 5, 4, 3, 2]
        #self.assertTrue(ogc_odd.perm_sign(G, p) == 1, 'incorrect permutation sign')
        #self.assertTrue(ogc_even.perm_sign(G, p) == -1, 'incorrect permutation sign')
        p = [0, 1, 2, 4, 3, 5]
        self.assertTrue(ogc_odd.perm_sign(G, p) == -1, 'incorrect permutation sign')
        self.assertTrue(ogc_even.perm_sign(G, p) == 1, 'incorrect permutation sign')

    def test_create_basis(self):

        for evenEdges in [True, False]:
            ogc = OGC.OrdinaryGraphVectorSpace(10, 9, evenEdges=evenEdges)
            self.GVSTest.test_create_basis(ogc)
            ogc_on_fly = OGC.OrdinaryGraphVectorSpace(6, 5, evenEdges=evenEdges, basis_on_fly=True)
            self.GVSTest.test_create_basis_on_fly(ogc,ogc_on_fly)

def suite():
    suite = unittest.TestSuite()
    suite.addTest(OGCTestCase('test_perm_sign'))
    suite.addTest(OGCTestCase('test_create_basis'))
    return suite

if __name__ == '__main__':
    runner = unittest.TextTestRunner()
    runner.run(suite())
