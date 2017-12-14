import unittest
import GraphVectorSpace as GVS
import OrdinaryGraphComplex as OGC

reload(OGC)


class GVSTestCase(unittest.TestCase):
    def setUp(self):
        pass

    def test_create_basis(self, graph_vector_space):

        graph_vector_space.delete_file()
        self.assertFalse(graph_vector_space.basis_built(), 'basis should not be built yet')
        basis1 = graph_vector_space.create_basis_g6()
        self.assertTrue(graph_vector_space.basis_built(), 'basis should be built')
        basis2 = graph_vector_space.get_basis(g6=True)
        self.assertListEqual(basis1, basis2, 'created and loaded basis are not equal')
        graph_vector_space.delete_file()
        self.assertFalse(graph_vector_space.basis_built(), 'basis file should have been deleted')

    def test_create_basis_on_fly(self, graph_vector_space, graph_vector_space_on_fly):

        self.assertFalse(graph_vector_space_on_fly.basis_built(), 'basis should be built')
        basis1 = graph_vector_space_on_fly.get_basis_g6()
        graph_vector_space_on_fly.delet_basis()
        self.assertFalse(graph_vector_space_on_fly.basis_built(), 'basis file should have been deleted')
        basis2 = graph_vector_space.get_basis()
        self.assertListEqual(basis1, basis2, 'created and loaded basis(created on fly) are not equal')
        graph_vector_space.delete_file()
        self.assertFalse(graph_vector_space.basis_built(), 'basis file should have been deleted')

