import unittest

from phylox import DiNetwork
from phylox.cherrypicking.combining_networks import HybridizationProblem
from phylox.cherrypicking.tree_child_sequences import check_cherry_picking_sequence
from phylox.constants import LENGTH_ATTR


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
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_trees_isomorphic(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "A"), (7, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_trees_different_but_trivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "C"), (7, "D")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        # resulting sequence should reduce both, but also be a valid sequence
        # hence, it also adds a root.
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_trees_different_nontrivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (4, 5), (4, 6)],
            labels=[(3, "A"), (5, "B"), (6, "C")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (4, 5), (4, 6)],
            labels=[(5, "A"), (3, "B"), (6, "C")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_three_trees(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "A"), (7, "B")],
        )
        network3 = DiNetwork(
            edges=[(1, 2), (2, 8), (2, 9)],
            labels=[(8, "A"), (9, "B")],
        )
        networks = [network1, network2, network3]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_one_network(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 2)
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_networks_trivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 7), (4, 8)],
            labels=[(7, "A"), (8, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 2)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_networks_nontrivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 7), (4, 8)],
            labels=[(8, "A"), (7, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True)
        print(result)
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_multiple_large_networks_same_labels(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 9),
                (6, 10),
            ],
            labels=[(7, "A"), (8, "B"), (9, "C"), (10, "D")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "A"), (7, "B"), (11, "C"), (12, "D")],
        )
        network3 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "A"), (11, "B"), (12, "C"), (7, "D")],
        )
        networks = [network1, network2, network3]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_multiple_large_networks_different_labels(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 9),
                (6, 10),
            ],
            labels=[(7, "A"), (8, "B"), (9, "C"), (10, "D")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "F"), (7, "B"), (11, "C"), (12, "E")],
        )
        network3 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "A"), (11, "E"), (12, "G"), (7, "F")],
        )
        networks = [network1, network2, network3]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(progress=True, seed=1)
        for network in networks:
            print(network.edges(data=True))
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_seed(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3,4), (3,5), (4,6), (5,6), (5,7), (6,8), (7,8), (7,9), (8,10), (9,10), (9,11), (10,12)],
            labels=[(11, "A"), (12, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3,4), (3,5), (4,6), (5,6), (5,7), (6,8), (7,8), (7,9), (8,10), (9,10), (9,11), (10,12)],
            labels=[(11, "B"), (12, "A")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result1 = problem.CPSBound(progress=True, seed=1)
        problem.best_seq = None
        problem.best_seq_with_lengths = None
        result2 = problem.CPSBound(progress=True, seed=1)
        problem.best_seq = None
        problem.best_seq_with_lengths = None
        result3 = problem.CPSBound(progress=True, seed=2)
        self.assertEqual(result1, result2)
        self.assertNotEqual(result1, result3)


class TestHybridizationProblemTrackPairs(unittest.TestCase):
    def test_one_tree(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_trees_isomorphic(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "A"), (7, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_trees_different_but_trivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "C"), (7, "D")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        # resulting sequence should reduce both, but also be a valid sequence
        # hence, it also adds a root.
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_trees_different_nontrivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (4, 5), (4, 6)],
            labels=[(3, "A"), (5, "B"), (6, "C")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (4, 5), (4, 6)],
            labels=[(5, "A"), (3, "B"), (6, "C")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_three_trees(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 6), (2, 7)],
            labels=[(6, "A"), (7, "B")],
        )
        network3 = DiNetwork(
            edges=[(1, 2), (2, 8), (2, 9)],
            labels=[(8, "A"), (9, "B")],
        )
        networks = [network1, network2, network3]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_one_network(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        self.assertEqual(len(result), 2)
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_networks_trivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 7), (4, 8)],
            labels=[(7, "A"), (8, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        self.assertEqual(len(result), 2)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_networks_nontrivial(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 7), (4, 8)],
            labels=[(8, "A"), (7, "B")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        print(result)
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_multiple_large_networks_same_labels(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 9),
                (6, 10),
            ],
            labels=[(7, "A"), (8, "B"), (9, "C"), (10, "D")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "A"), (7, "B"), (11, "C"), (12, "D")],
        )
        network3 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "A"), (11, "B"), (12, "C"), (7, "D")],
        )
        networks = [network1, network2, network3]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_multiple_large_networks_different_labels(self):
        network1 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 9),
                (6, 10),
            ],
            labels=[(7, "A"), (8, "B"), (9, "C"), (10, "D")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "F"), (7, "B"), (11, "C"), (12, "E")],
        )
        network3 = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 7),
                (5, 8),
                (6, 8),
                (6, 10),
                (8, 9),
                (9, 11),
                (9, 12),
            ],
            labels=[(10, "A"), (11, "E"), (12, "G"), (7, "F")],
        )
        networks = [network1, network2, network3]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True, seed=1)
        for network in networks:
            print(network.edges(data=True))
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))
        
    def test_seed(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3,4), (3,5), (4,6), (5,6), (5,7), (6,8), (7,8), (7,9), (8,10), (9,10), (9,11), (10,12)],
            labels=[(11, "A"), (12, "B")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3,4), (3,5), (4,6), (5,6), (5,7), (6,8), (7,8), (7,9), (8,10), (9,10), (9,11), (10,12)],
            labels=[(11, "B"), (12, "A")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result1 = problem.CPSBound(track=True, progress=True, seed=1)
        problem.best_seq = None
        problem.best_seq_with_lengths = None
        result2 = problem.CPSBound(track=True, progress=True, seed=1)
        problem.best_seq = None
        problem.best_seq_with_lengths = None
        result3 = problem.CPSBound(track=True, progress=True, seed=2)
        print(result1)
        print(result2)
        print(result3)
        self.assertEqual(result1, result2)
        self.assertNotEqual(result1, result3)


class TestHybridizationProblemLengths(unittest.TestCase):
    def test_one_tree(self):
        network = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(3, "A"), (4, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        self.assertEqual(len(result), 1)
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_multiple_missing_lengths(self):
        network1 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(3, "A"), (4, "B")],
        )
        network3 = DiNetwork(
            edges=[(1, 2), (2, 8), (2, 9)],
            labels=[(8, "A"), (9, "B")],
        )
        with self.assertRaises(ValueError):
            problem = HybridizationProblem(
                [network1, network2, network3], newick_strings=False
            )
            result = problem.CPSBound(lengths=True, progress=True)

    def test_two_trees_isomorphic(self):
        network1 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 1.5}),
                (2, 3, {LENGTH_ATTR: 1.2}),
                (2, 4, {LENGTH_ATTR: 0.8}),
            ],
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
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(3, "A"), (4, "B")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
            ],
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
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
                (3, 4, {LENGTH_ATTR: 1.0}),
                (3, 5, {LENGTH_ATTR: 1.0}),
                (4, 6, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(5, "A"), (6, "B")],
        )
        problem = HybridizationProblem([network], newick_strings=False)
        result = problem.CPSBound(lengths=True, progress=True)
        print(result)
        self.assertEqual(len(result), 2)
        self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))

    def test_two_networks_trivial(self):
        network1 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
                (3, 4, {LENGTH_ATTR: 1.0}),
                (3, 5, {LENGTH_ATTR: 1.0}),
                (4, 6, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 0.8}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.4}),
                (3, 4, {LENGTH_ATTR: 1.1}),
                (3, 5, {LENGTH_ATTR: 1.0}),
                (4, 6, {LENGTH_ATTR: 1.6}),
            ],
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
            edges=[
                (1, 2, {LENGTH_ATTR: 1.0}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.0}),
                (3, 4, {LENGTH_ATTR: 1.0}),
                (3, 5, {LENGTH_ATTR: 1.0}),
                (4, 6, {LENGTH_ATTR: 1.0}),
            ],
            labels=[(5, "A"), (6, "B")],
        )
        network2 = DiNetwork(
            edges=[
                (1, 2, {LENGTH_ATTR: 0.8}),
                (2, 3, {LENGTH_ATTR: 1.0}),
                (2, 4, {LENGTH_ATTR: 1.4}),
                (3, 4, {LENGTH_ATTR: 1.1}),
                (3, 5, {LENGTH_ATTR: 1.0}),
                (4, 6, {LENGTH_ATTR: 1.6}),
            ],
            labels=[(5, "B"), (6, "A")],
        )
        networks = [network1, network2]
        problem = HybridizationProblem(networks, newick_strings=False)
        result = problem.CPSBound(track=True, progress=True, seed=1)
        self.assertEqual(len(result), 3)
        for network in networks:
            self.assertTrue(check_cherry_picking_sequence(network, result, labels=True))