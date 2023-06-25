import unittest

from phylox import DiNetwork
from phylox.cherrypicking.tree_child_sequences import *


class TestTreeChildContainment(unittest.TestCase):
    def test_simple_cherry_unlabeled(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
        )
        result = tree_child_network_contains(network, network)
        self.assertTrue(result)
    
    def test_simple_cherry_picking_unlabeled(self):
        network1 = DiNetwork(
            edges=[(0, 1), (1, 2), (1, 3), (2, 3), (2,4), (3, 5)],
        )
        network2 = DiNetwork(
            edges=[(0, 1), (1, 4), (1, 5)],
        )
        result = tree_child_network_contains(network1, network2)
        self.assertTrue(result)
        #repeat to check that networks are preserved
        result = tree_child_network_contains(network1, network2)
        self.assertTrue(result)
        result = tree_child_network_contains(network2, network1)
        self.assertFalse(result)
