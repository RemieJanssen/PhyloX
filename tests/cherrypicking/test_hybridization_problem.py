import unittest

from phylox import DiNetwork
from phylox.cherrypicking.combining_networks import HybridizationProblem


class TestHybridizationProblem(unittest.TestCase):
    def test_one_tree(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound()
        print(result)
        self.assertEqual(len(result), 1)