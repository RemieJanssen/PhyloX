import unittest

from phylox.constants import LABEL_ATTR
from phylox.generators.randomTC import (
    generate_network_random_tree_child_sequence,
    random_tree_child_sequence,
    random_tree_child_subsequence,
)
from phylox import DiNetwork


class TestRandomTC(unittest.TestCase):
    def test_one_leaf(self):
        with self.assertRaises(ValueError):
            generate_network_random_tree_child_sequence(
                leaves=1,
                reticulations=0,
            )

    def test_two_leaves(self):
        leaves = 2
        reticulations = 0
        network = generate_network_random_tree_child_sequence(
            leaves=leaves,
            reticulations=reticulations,
        )
        self.assertEqual(len(network.leaves), leaves)
        self.assertEqual(len(network.roots), 1)
        self.assertEqual(network.reticulation_number, reticulations)

    def test_tree(self):
        leaves = 5
        reticulations = 0
        network = generate_network_random_tree_child_sequence(
            leaves=leaves,
            reticulations=reticulations,
        )
        self.assertEqual(len(network.leaves), leaves)
        self.assertEqual(len(network.roots), 1)
        self.assertEqual(network.reticulation_number, reticulations)

    def test_network(self):
        leaves = 10
        reticulations = 5
        network = generate_network_random_tree_child_sequence(
            leaves=leaves,
            reticulations=reticulations,
        )
        self.assertEqual(len(network.leaves), leaves)
        self.assertEqual(len(network.roots), 1)
        self.assertEqual(network.reticulation_number, reticulations)

    def test_network_by_label(self):
        leaves = 10
        reticulations = 5
        network = generate_network_random_tree_child_sequence(
            leaves=leaves,
            reticulations=reticulations,
            label_leaves=True,
        )
        self.assertEqual(len(network.leaves), leaves)
        self.assertEqual(len(network.roots), 1)
        self.assertEqual(network.reticulation_number, reticulations)
        for node in network.leaves:
            self.assertIsNotNone(network.nodes[node].get(LABEL_ATTR))

    def test_subsequence(self):
        leaves = 10
        reticulations = 5
        sequence = random_tree_child_sequence(
            leaves=leaves,
            reticulations=reticulations,
        )
        network = DiNetwork.from_cherry_picking_sequence(sequence)
        self.assertEqual(len(network.leaves), leaves)
        self.assertEqual(len(network.roots), 1)
        self.assertEqual(network.reticulation_number, reticulations)
        subsequence_reticulations = 3
        subsequence = random_tree_child_subsequence(sequence, subsequence_reticulations)
        subnetwork = DiNetwork.from_cherry_picking_sequence(subsequence)
        self.assertEqual(len(subnetwork.leaves), leaves)
        self.assertEqual(len(subnetwork.roots), 1)
        self.assertEqual(subnetwork.reticulation_number, subsequence_reticulations)
