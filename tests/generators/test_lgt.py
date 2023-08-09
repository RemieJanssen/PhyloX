import unittest

from phylox.generators.lgt import generate_network_lgt
from phylox.isomorphism import is_isomorphic

class TestLGTNetwork(unittest.TestCase):
    def test_one_leaf(self):
        network = generate_network_lgt(
            n=1,
            k=0,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 1)
        self.assertEqual(len(network.edges), 1)

    def test_one_leaf_many_retics(self):
        network = generate_network_lgt(
            n=1,
            k=100,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 1)
        self.assertEqual(len(network.reticulations), 100)

    def test_two_leaves(self):
        network = generate_network_lgt(
            n=2,
            k=0,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 2)
        self.assertEqual(len(network.edges), 3)

    def test_one_reticulation(self):
        network = generate_network_lgt(
            n=4,
            k=1,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 4)
        self.assertEqual(len(network.reticulations), 1)

    def test_larger_internal(self):
        network = generate_network_lgt(
            n=100,
            k=10,
            wint=1,
            wext=0.001,
        )
        self.assertEqual(len(network.leaves), 100)
        self.assertEqual(len(network.reticulations), 10)

    def test_larger_external(self):
        network = generate_network_lgt(
            n=100,
            k=10,
            wint=0.001,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 100)
        self.assertEqual(len(network.reticulations), 10)

    def test_seed(self):
        network1 = generate_network_lgt(
            n=100,
            k=10,
            wint=1,
            wext=1,
            seed=1,
        )
        network2 = generate_network_lgt(
            n=100,
            k=10,
            wint=1,
            wext=1,
            seed=1,
        )
        network3 = generate_network_lgt(
            n=100,
            k=10,
            wint=1,
            wext=1,
            seed=2,
        )
        self.assertTrue(is_isomorphic(network1, network2))
        self.assertFalse(is_isomorphic(network1, network3))
