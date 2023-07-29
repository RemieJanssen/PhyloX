import unittest

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR
from phylox.parser import extended_newick_to_dinetwork, dinetwork_to_extended_newick
from phylox.isomorphism import is_isomorphic


class TestExtendedNewickToDiNetwork(unittest.TestCase):
    def test_small_tree(self):
        newick = "(a,b,(c,d));"
        network = extended_newick_to_dinetwork(newick)
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (1, 4), (4, 5), (4, 6)],
            labels=[(2, "a"), (3, "b"), (5, "c"), (6, "d")],
        )
        self.assertTrue(is_isomorphic(network, network2))

    def test_small_tree_with_lengths(self):
        newick = "(a:1.0,b:1.1,(c:1.2,d:1.3):1.4);"
        network = extended_newick_to_dinetwork(newick)
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (1, 4), (4, 5), (4, 6)],
            labels=[(2, "a"), (3, "b"), (5, "c"), (6, "d")],
        )
        self.assertTrue(is_isomorphic(network, network2))
        node_a = network.label_to_node_dict["a"]
        parent_a = network.parent(node_a)
        self.assertEqual(network[parent_a][node_a]["length"], 1.0)

    def test_small_tree_with_more_attrs(self):
        newick = "(a:1.0:2.0:3.0,b:1.1,(c:1.2,d:1.3):1.4);"
        network = extended_newick_to_dinetwork(newick)
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (1, 4), (4, 5), (4, 6)],
            labels=[(2, "a"), (3, "b"), (5, "c"), (6, "d")],
        )
        self.assertTrue(is_isomorphic(network, network2))
        node_a = network.label_to_node_dict["a"]
        parent_a = network.parent(node_a)
        self.assertEqual(network[parent_a][node_a]["length"], 1.0)
        self.assertEqual(network[parent_a][node_a]["bootstrap"], 2.0)
        self.assertEqual(network[parent_a][node_a]["probability"], 3.0)

    def test_network(self):
        newick = "(a,(b)#R1,(#H1,c));"
        network = extended_newick_to_dinetwork(newick)
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (3, 7), (1, 4), (4, 3), (4, 6)],
            labels=[(2, "a"), (7, "b"), (6, "c")],
        )
        self.assertTrue(is_isomorphic(network, network2))

    def test_network_with_lengths(self):
        newick = "(a:1.0,(b:1.1)#R1:1.2,(#H1:1.3,c:1.4):1.5);"
        network = extended_newick_to_dinetwork(newick)
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (3, 7), (1, 4), (4, 3), (4, 6)],
            labels=[(2, "a"), (7, "b"), (6, "c")],
        )
        self.assertTrue(is_isomorphic(network, network2))
        node_a = network.label_to_node_dict["a"]
        parent_a = network.parent(node_a)
        self.assertEqual(network[parent_a][node_a]["length"], 1.0)


class TestNetworkToNewick(unittest.TestCase):
    def test_small_tree(self):
        network = DiNetwork(
            edges=[(1, 2), (1, 3), (1, 4), (4, 5), (4, 6)],
            labels=[(2, "a"), (3, "b"), (5, "c"), (6, "d")],
        )
        newick = dinetwork_to_extended_newick(network)
        assert False
