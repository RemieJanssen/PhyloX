import unittest

# from phylox import DiNetwork
from phylox.generators.trees.beta_splitting_tree import simulate_beta_splitting


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

