import unittest

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR


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
