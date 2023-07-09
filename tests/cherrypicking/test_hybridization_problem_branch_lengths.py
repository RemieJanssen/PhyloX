import unittest

from phylox import DiNetwork
from phylox.cherrypicking.combining_networks import HybridizationProblem
from phylox.cherrypicking.tree_child_sequences import check_cherry_picking_sequence
from phylox.constants import LENGTH_ATTR

class TestHybridizationProblem(unittest.TestCase):
    def test_one_tree(self):
        network = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0})],
            labels=[(3, "A"), (4, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_multiple_missing_lengths(self):
        network1 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0})],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0})],
            labels=[(3, "A"), (4, "B")],
        )
        network3 = DiNetwork(
            edges=[(1, 2), (2, 8), (2, 9)],
            labels=[(8, "A"), (9, "B")],
        )
        with self.assertRaises(ValueError):
            problem = HybridizationProblem([network1, network2, network3], newick_strings=False)
            result = problem.CPSBound(lengths=True, progress=True)

    def test_two_trees_isomorphic(self):
        network1 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0})],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.5}), (2, 3, {LENGTH_ATTR: 1.2}), (2, 4, {LENGTH_ATTR: 0.8})],
            labels=[(3, "A"), (4, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_trees_different_but_trivial(self):
        network1 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0})],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0})],
            labels=[(3, "C"), (4, "D")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        # resulting sequence should reduce both, but also be a valid sequence
        # hence, it also adds a root.
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_one_network(self):
        network = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0}), (3,4, {LENGTH_ATTR: 1.0}), (3,5, {LENGTH_ATTR: 1.0}), (4,6, {LENGTH_ATTR: 1.0})],
            labels=[(5, "A"), (6, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        self.assertEqual(len(result), 2)
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_networks_trivial(self):
        network1 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0}), (3,4, {LENGTH_ATTR: 1.0}), (3,5, {LENGTH_ATTR: 1.0}), (4,6, {LENGTH_ATTR: 1.0})],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 0.8}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.4}), (3,4, {LENGTH_ATTR: 1.1}), (3,5, {LENGTH_ATTR: 1.0}), (4,6, {LENGTH_ATTR: 1.6})],
            labels=[(5, "A"), (6, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        self.assertEqual(len(result), 2)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))


    def test_two_networks_nontrivial(self):
        network1 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 1.0}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.0}), (3,4, {LENGTH_ATTR: 1.0}), (3,5, {LENGTH_ATTR: 1.0}), (4,6, {LENGTH_ATTR: 1.0})],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2, {LENGTH_ATTR: 0.8}), (2, 3, {LENGTH_ATTR: 1.0}), (2, 4, {LENGTH_ATTR: 1.4}), (3,4, {LENGTH_ATTR: 1.1}), (3,5, {LENGTH_ATTR: 1.0}), (4,6, {LENGTH_ATTR: 1.6})],
            labels=[(5, "B"), (6, "A")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))
