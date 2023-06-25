import unittest

from phylox import DiNetwork
from phylox.networkproperties.properties import *


class TestCountReduciblePairs(unittest.TestCase):
    def test_simple_cherry(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        result = count_reducible_pairs(network)
        self.assertEqual(result["cherries"], 1)
        self.assertEqual(result["reticulate_cherries"], 0)

    def test_no_network(self):
        network = DiNetwork()
        result = count_reducible_pairs(network)
        self.assertEqual(result["cherries"], 0)
        self.assertEqual(result["reticulate_cherries"], 0)

    def test_simple_reticulated_cherry(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        result = count_reducible_pairs(network)
        self.assertEqual(result["cherries"], 0)
        self.assertEqual(result["reticulate_cherries"], 1)


class TestBlobProperties(unittest.TestCase):
    def test_no_blob(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        result = blob_properties(network)
        self.assertEqual(result, [])

    def test_simple_blob(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        result = blob_properties(network)
        self.assertEqual(result, [(3, 1)])

    def test_one_large_blob(self):
        network = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (5, 6),
                (5, 7),
                (6, 8),
                (7, 8),
                (7, 9),
                (8, 10),
            ],
        )
        result = blob_properties(network)
        self.assertEqual(result, [(7, 3)])

    def test_two_blobs(self):
        network = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 4),
                (3, 5),
                (4, 6),
                (6, 7),
                (6, 8),
                (7, 8),
                (7, 9),
                (8, 10),
            ],
        )
        result = blob_properties(network)
        self.assertEqual(result, [(3, 1), (3, 1)])


class TestB2Balance(unittest.TestCase):
    def test_no_network(self):
        network = DiNetwork()
        result = b2_balance(network)
        self.assertEqual(result, 0)

    def test_path(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (3, 4)],
        )
        result = b2_balance(network)
        self.assertEqual(result, 0)

    def test_simple_network(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        result = b2_balance(network)
        self.assertEqual(result, 1)

    def test_balanced_four_leaves(self):
        network = DiNetwork(
            edges=[(1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7)],
        )
        result = b2_balance(network)
        self.assertEqual(result, 2)

    def test_balanced_four_leaves_with_root_edges(self):
        network = DiNetwork(
            edges=[(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7)],
        )
        result = b2_balance(network)
        self.assertEqual(result, 2)

    def test_caterpillar_four_leaves(self):
        network = DiNetwork(
            edges=[(0, 1), (1, 2), (1, 3), (3, 4), (3, 5), (5, 6), (5, 7)],
        )
        result = b2_balance(network)
        self.assertEqual(result, 0.5 + 2 / 4 + 6 / 8)

    def test_connect_roots(self):
        network = DiNetwork(
            edges=[(0, 1), (1, 2), (1, 3), (4, 5), (5, 6), (5, 7)],
        )
        with self.assertRaises(ValueError):
            b2_balance(network, connect_roots=False)
        result = b2_balance(network, connect_roots=True)
        self.assertEqual(result, 2)
