import unittest

from phylox.generators.zods import generate_network_zods


class TestZodsNetwork(unittest.TestCase):
    def test_no_time(self):
        network = generate_network_zods(
            time_limit=0,
            speciation_rate=1,
            hybridization_rate=0,
        )
        self.assertEqual(len(network.leaves), 1)
        self.assertEqual(len(network.edges), 1)
        self.assertEqual(len(network.roots), 1)

    def test_no_reticulations(self):
        network = generate_network_zods(
            time_limit=3,
            speciation_rate=1,
            hybridization_rate=0,
        )
        leaves = len(network.leaves)
        self.assertEqual(len(network.reticulations), 0)
        self.assertGreater(leaves, 0)
        self.assertEqual(len(network.edges), 2 * leaves - 1)
        self.assertEqual(len(network.roots), 1)

    def test_some_reticulations(self):
        network = generate_network_zods(
            time_limit=50,
            speciation_rate=1,
            hybridization_rate=0.5,
        )
        print(network.edges)
        leaves = len(network.leaves)
        reticulations = len(network.reticulations)
        self.assertGreater(reticulations, 0)
        self.assertGreater(leaves, 0)
        self.assertEqual(len(network.edges), 2 * leaves + 3 * reticulations - 1)
        self.assertEqual(len(network.roots), 1)
