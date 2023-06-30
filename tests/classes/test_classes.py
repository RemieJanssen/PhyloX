import unittest

from phylox import DiNetwork
from phylox.classes import *
from phylox.constants import LABEL_ATTR


class TestClassBinary(unittest.TestCase):
    def test_is_binary(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        self.assertTrue(is_binary(network1))

    def test_is_binary_false_tree_node(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (2, 5)],
        )
        self.assertFalse(is_binary(network1))

    def test_is_binary_false_retic(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (4, 5),
                (4, 6),
                (3, 7),
                (5, 7),
                (6, 7),
                (6, 8),
                (7, 9),
            ],
        )
        self.assertFalse(is_binary(network1))

    def test_is_binary_false_mixed(self):
        # 4 is a mixed node with indegree 2 and outdegree 2
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6), (4, 7)],
        )
        self.assertFalse(is_binary(network1))


class TestClassTreeChild(unittest.TestCase):
    def test_is_tree_child(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        self.assertTrue(is_tree_child(network1))

    def test_is_tree_child_2(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        self.assertTrue(is_tree_child(network1))

    def test_is_tree_child_false_stack(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (4, 5),
                (3, 8),
                (8, 5),
                (8, 9),
                (5, 6),
            ],
        )
        self.assertFalse(is_tree_child(network1))

    def test_is_tree_child_false_w_shape(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 8),
                (8, 5),
                (8, 9),
                (5, 6),
            ],
        )
        self.assertFalse(is_tree_child(network1))

    def test_is_tree_child_false_diamond(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (4, 5), (3, 5), (5, 6)],
        )
        self.assertFalse(is_tree_child(network1))


class TestClassOrchard(unittest.TestCase):
    def test_is_orchard_empty_network(self):
        network1 = DiNetwork()
        self.assertTrue(is_orchard(network1))

    def test_is_orchard(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        self.assertTrue(is_orchard(network1))

    def test_is_orchard_add_root(self):
        network1 = DiNetwork(
            edges=[(2, 3), (2, 4), (3, 4), (3, 5), (4, -1)],
        )
        self.assertTrue(is_orchard(network1))

    def test_is_orchard_false(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 5),
                (4, 5),
                (3, 6),
                (4, 6),
                (5, 7),
                (6, 8),
            ],
        )
        self.assertFalse(is_orchard(network1))


class TestClassStackFree(unittest.TestCase):
    def test_is_stack_free(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        self.assertTrue(is_stack_free(network1))

    def test_is_stack_free_2(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        self.assertTrue(is_stack_free(network1))

    def test_is_stack_free_false_stack(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (4, 5),
                (3, 8),
                (8, 5),
                (8, 9),
                (5, 6),
            ],
        )
        self.assertFalse(is_stack_free(network1))

    def test_is_stack_free_w_shape(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 8),
                (8, 5),
                (8, 9),
                (5, 6),
            ],
        )
        self.assertTrue(is_stack_free(network1))

    def test_is_stack_free_false_diamond(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (4, 5), (3, 5), (5, 6)],
        )
        self.assertFalse(is_stack_free(network1))


class TestClassTreeBased(unittest.TestCase):
    def test_is_tree_based(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        self.assertTrue(is_tree_based(network1))

    def test_is_tree_based_2(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        self.assertTrue(is_tree_based(network1))

    def test_is_tree_based_stack(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (4, 5),
                (3, 8),
                (8, 5),
                (8, 9),
                (5, 6),
            ],
        )
        self.assertTrue(is_tree_based(network1))

    def test_is_tree_based_w_shape(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 8),
                (8, 5),
                (8, 9),
                (5, 6),
            ],
        )
        self.assertTrue(is_tree_based(network1))

    def test_is_tree_based_diamond(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (4, 5), (3, 5), (5, 6)],
        )
        self.assertTrue(is_tree_based(network1))

    def test_is_tree_based_false(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (3, 4),
                (3, 5),
                (4, 5),
                (4, 6),
                (5, 7),
                (6, 7),
                (7, 8),
                (2, 6),
            ],
        )
        self.assertFalse(is_tree_based(network1))
