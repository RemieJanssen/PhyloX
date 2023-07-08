import unittest

from phylox import DiNetwork
from phylox.cherrypicking.combining_networks import HybridizationProblem


class TestHybridizationProblem(unittest.TestCase):
    def test_one_tree(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 1)

    def test_two_trees_isomorphic(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "A"), (7, "B")],
        )
        problem = HybridizationProblem([network1, network2], newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 1)

    def test_two_trees_different_but_trivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "C"), (7, "D")],
        )
        problem = HybridizationProblem([network1, network2], newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 2)

    def test_two_trees_different_nontrivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (4, 5), (4, 6)],
            labels=[(3, "A"), (5, "B"), (6, "C")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (4, 5), (4, 6)],
            labels=[(5, "A"), (3, "B"), (6, "C")],
        )
        problem = HybridizationProblem([network1, network2], newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 3)