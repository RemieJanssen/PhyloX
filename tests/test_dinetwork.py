import unittest

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR
from phylox.isomorphism import is_isomorphic


class TestDiNetwork(unittest.TestCase):
    def test_init(self):
        network = DiNetwork()
        self.assertEqual(list(network.edges), [])
        self.assertEqual(list(network.nodes), [])

    def test_init_all(self):
        network = DiNetwork(
            nodes=[1, 2, 3],
            edges=[(1, 2), (2, 3)],
            labels=[(1, "a"), (2, "b"), (3, "c")],
        )
        self.assertCountEqual(list(network.nodes), [1, 2, 3])
        self.assertCountEqual(list(network.edges), [(1, 2), (2, 3)])
        self.assertEqual(network.nodes[1][LABEL_ATTR], "a")
        self.assertEqual(network.nodes[2][LABEL_ATTR], "b")
        self.assertEqual(network.nodes[3][LABEL_ATTR], "c")

    def test_leaves(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3)],
            labels=[(1, "a"), (2, "b"), (3, "c")],
        )
        self.assertEqual(network.leaves, {3})

    def test_roots(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3)],
            labels=[(1, "a"), (2, "b"), (3, "c")],
        )
        self.assertEqual(network.roots, {1})

    def test_reticulation_number(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertEqual(network.reticulation_number, 1)

    def test_child(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertIn(network.child(2), [3, 4])
        self.assertEqual(network.child(2, exclude=[3]), 4)
        self.assertEqual(network.child(2, exclude=[3, 4]), None)
        self.assertIn(network.child(2, randomNodes=True), [3, 4])
        self.assertEqual(network.child(2, exclude=[3], randomNodes=True), 4)
        self.assertEqual(network.child(2, exclude=[3, 4], randomNodes=True), None)

    def test_child_seed(self):
        network = DiNetwork(
            edges=[(0, i) for i in range(1, 100)],
        )
        child1 = network.child(0, randomNodes=True, seed=1)
        child2 = network.child(0, randomNodes=True, seed=1)
        child3 = network.child(0, randomNodes=True, seed=2)
        self.assertEqual(child1, child2)
        self.assertNotEqual(child1, child3)

    def test_parent(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertIn(network.parent(4), [2, 3])
        self.assertEqual(network.parent(4, exclude=[2]), 3)
        self.assertEqual(network.parent(4, exclude=[2, 3]), None)
        self.assertIn(network.parent(4, randomNodes=True), [2, 3])
        self.assertEqual(network.parent(4, exclude=[2], randomNodes=True), 3)
        self.assertEqual(network.parent(4, exclude=[2, 3], randomNodes=True), None)

    def test_parent_seed(self):
        network = DiNetwork(
            edges=[(i, 0) for i in range(1, 100)],
        )
        parent1 = network.parent(0, randomNodes=True, seed=1)
        parent2 = network.parent(0, randomNodes=True, seed=1)
        parent3 = network.parent(0, randomNodes=True, seed=2)
        self.assertEqual(parent1, parent2)
        self.assertNotEqual(parent1, parent3)

    def test_is_tree_node(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertFalse(network.is_tree_node(1))
        self.assertTrue(network.is_tree_node(2))
        self.assertTrue(network.is_tree_node(3))
        self.assertFalse(network.is_tree_node(4))
        self.assertFalse(network.is_tree_node(5))
        self.assertFalse(network.is_tree_node(6))

    def test_is_reticulation(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertFalse(network.is_reticulation(1))
        self.assertFalse(network.is_reticulation(2))
        self.assertFalse(network.is_reticulation(3))
        self.assertTrue(network.is_reticulation(4))
        self.assertFalse(network.is_reticulation(5))
        self.assertFalse(network.is_reticulation(6))

    def test_is_leaf(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertFalse(network.is_leaf(1))
        self.assertFalse(network.is_leaf(2))
        self.assertFalse(network.is_leaf(3))
        self.assertFalse(network.is_leaf(4))
        self.assertTrue(network.is_leaf(5))
        self.assertTrue(network.is_leaf(6))

    def test_is_root(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertTrue(network.is_root(1))
        self.assertFalse(network.is_root(2))
        self.assertFalse(network.is_root(3))
        self.assertFalse(network.is_root(4))
        self.assertFalse(network.is_root(5))
        self.assertFalse(network.is_root(6))

    def test_from_newick(self):
        network = DiNetwork.from_newick("((a,b),c);")
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (2, 4), (2, 5)],
            labels=[(3, "c"), (4, "b"), (5, "a")],
        )
        self.assertTrue(is_isomorphic(network, network2))

    def test_to_newick(self):
        network = DiNetwork.from_newick("(a,b);")
        newick = network.newick()
        self.assertTrue(newick in ["(a,b);", "(b,a);"])
