import unittest

# from phylox import DiNetwork
from phylox.generators.trees.beta_splitting_tree import simulate_beta_splitting
from phylox.generators.trees.well_known import generate_balanced_tree, generate_caterpillar
from phylox.generators.trees.add_edges import network_from_tree, AddEdgeMethod


class TestBetaSplitting(unittest.TestCase):
    def test_one_leaf(self):
        tree = simulate_beta_splitting(
            n=1,
            beta=1,
        )
        self.assertEqual(len(tree.leaves), 1)
        self.assertEqual(len(tree.edges), 1)

    def test_two_leaves(self):
        tree = simulate_beta_splitting(
            n=2,
            beta=.5,
        )
        self.assertEqual(len(tree.leaves), 2)
        self.assertEqual(len(tree.edges), 3)
    
    def test_larger(self):
        tree = simulate_beta_splitting(
            n=10,
            beta=0,
        )
        self.assertEqual(len(tree.leaves), 10)
        self.assertEqual(len(tree.edges), 2*10-1)

class TestWellKnown(unittest.TestCase):
    def test_balanced_tree(self):
        tree = generate_balanced_tree(8)
        self.assertEqual(len(tree.leaves), 8)
        self.assertEqual(len(tree.edges), 2*8-1)

    def test_caterpillar(self):
        number_of_leaves = 5
        tree = generate_caterpillar(number_of_leaves)
        self.assertEqual(len(tree.leaves), number_of_leaves)
        self.assertEqual(len(tree.edges), 2*number_of_leaves-1)

class TestAddEdges(unittest.TestCase):
    def test_add_edges_bottom(self):
        tree = generate_balanced_tree(8)
        network = network_from_tree(tree, 3, AddEdgeMethod.BOTTOM)
        self.assertEqual(len(network.leaves), 8)
        self.assertEqual(len(network.edges), 2*8+3*3-1)

    def test_add_edges_uniform(self):
        tree = generate_balanced_tree(8)
        network = network_from_tree(tree, 3, AddEdgeMethod.UNIFORM)
        self.assertEqual(len(network.leaves), 8)
        self.assertEqual(len(network.edges), 2*8+3*3-1)
    
    def test_add_edges_local(self):
        tree = generate_balanced_tree(8)
        network = network_from_tree(tree, 3, AddEdgeMethod.LOCAL)
        self.assertEqual(len(network.leaves), 8)
        self.assertEqual(len(network.edges), 2*8+3*3-1)

    def test_add_edges_non_existent(self):
        tree = generate_balanced_tree(8)
        with self.assertRaises(ValueError):
            network = network_from_tree(tree, 3, "TOP")
