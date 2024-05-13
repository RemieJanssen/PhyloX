import unittest

from phylox.generators.lgt import generate_network_lgt, generate_network_lgt_conditional
from phylox.isomorphism import is_isomorphic


class TestLGTNetwork(unittest.TestCase):
    def test_one_leaf(self):
        network = generate_network_lgt(
            leaves_goal=1,
            retics_goal=0,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 1)
        self.assertEqual(len(network.edges), 1)

    def test_one_leaf_many_retics(self):
        network = generate_network_lgt(
            leaves_goal=1,
            retics_goal=100,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 1)
        self.assertEqual(len(network.reticulations), 100)

    def test_two_leaves(self):
        network = generate_network_lgt(
            leaves_goal=2,
            retics_goal=0,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 2)
        self.assertEqual(len(network.edges), 3)

    def test_one_reticulation(self):
        network = generate_network_lgt(
            leaves_goal=4,
            retics_goal=1,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 4)
        self.assertEqual(len(network.reticulations), 1)

    def test_larger_internal(self):
        network = generate_network_lgt(
            leaves_goal=100,
            retics_goal=10,
            wint=1,
            wext=0.001,
        )
        self.assertEqual(len(network.leaves), 100)
        self.assertEqual(len(network.reticulations), 10)

    def test_larger_external(self):
        network = generate_network_lgt(
            leaves_goal=100,
            retics_goal=10,
            wint=0.001,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 100)
        self.assertEqual(len(network.reticulations), 10)

    def test_seed(self):
        network1 = generate_network_lgt(
            leaves_goal=100,
            retics_goal=10,
            wint=1,
            wext=1,
            seed=1,
        )
        network2 = generate_network_lgt(
            leaves_goal=100,
            retics_goal=10,
            wint=1,
            wext=1,
            seed=1,
        )
        network3 = generate_network_lgt(
            leaves_goal=100,
            retics_goal=10,
            wint=1,
            wext=1,
            seed=2,
        )
        self.assertTrue(is_isomorphic(network1, network2))
        self.assertFalse(is_isomorphic(network1, network3))


class TestLGTNetworkConditional(unittest.TestCase):
    def test_one_leaf(self):
        network = generate_network_lgt_conditional(
            n=1,
            k=0,
            prob_lgt=0.5,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 1)
        self.assertEqual(len(network.edges), 1)

    def test_two_leaves(self):
        network = generate_network_lgt_conditional(
            n=2,
            k=0,
            prob_lgt=0.5,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 2)
        self.assertEqual(len(network.edges), 3)

    def test_one_reticulation(self):
        network = generate_network_lgt_conditional(
            n=4,
            k=1,
            prob_lgt=0.2,
            wint=1,
            wext=1,
        )
        self.assertEqual(len(network.leaves), 4)
        self.assertEqual(len(network.reticulations), 1)

    def test_seed(self):
        network1 = generate_network_lgt_conditional(
            n=4,
            k=1,
            prob_lgt=0.2,
            wint=1,
            wext=1,
            seed=1,
        )
        network2 = generate_network_lgt_conditional(
            n=4,
            k=1,
            prob_lgt=0.2,
            wint=1,
            wext=1,
            seed=1,
        )
        network3 = generate_network_lgt_conditional(
            n=4,
            k=1,
            prob_lgt=0.2,
            wint=1,
            wext=1,
            seed=2,
        )
        self.assertTrue(is_isomorphic(network1, network2))
        self.assertFalse(is_isomorphic(network1, network3))
