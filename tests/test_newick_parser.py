import unittest

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR
from phylox.isomorphism import is_isomorphic
from phylox.newick_parser import dinetwork_to_extended_newick, extended_newick_to_dinetwork


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
        newick = "(a:::3.0,b:1.1,(c,d:1.3)::1.4:);"
        network = extended_newick_to_dinetwork(newick)
        print(network.edges)
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (1, 4), (4, 5), (4, 6)],
            labels=[(2, "a"), (3, "b"), (5, "c"),(6, "d")],
        )
        self.assertTrue(is_isomorphic(network, network2))
        node_a = network.label_to_node_dict["a"]
        parent_a = network.parent(node_a)
        self.assertEqual(network[parent_a][node_a].get("length"), None)
        self.assertEqual(network[parent_a][node_a].get("bootstrap"), None)
        self.assertEqual(network[parent_a][node_a]["probability"], 3.0)

    def test_small_tree_with_partial_attrs(self):
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
    def test_multirooted(self):
        network = DiNetwork(
            edges=[(0, 1), (2, 3), (3, 4), (3, 5)],
            labels=[(1, "a"), (4, "b"), (5, "c")],
        )
        with self.assertRaises(ValueError):
            dinetwork_to_extended_newick(network)

    def test_small_tree(self):
        network = DiNetwork(
            edges=[(1, 2), (1, 3)],
            labels=[(2, "a"), (3, "b")],
        )
        newick = dinetwork_to_extended_newick(network)
        self.assertTrue(
            newick
            in [
                "(a,b);",
                "(b,a);",
            ]
        )

    def test_small_tree_lengths(self):
        network = DiNetwork(
            edges=[
                (1, 2, {"length": 1.0}),
                (1, 3, {"length": 2.0}),
            ],
            labels=[(2, "a"), (3, "b")],
        )
        newick = dinetwork_to_extended_newick(network)
        self.assertTrue(
            newick
            in [
                "(a:1.0,b:2.0);",
                "(b:2.0,a:1.0);",
            ]
        )

    def test_small_tree_more_attrs(self):
        network = DiNetwork(
            edges=[
                (1, 2, {"length": 1.0, "bootstrap": 2.0, "probability": 3.0}),
                (1, 3, {"length": 2.0, "bootstrap": 3.0, "probability": 4.0}),
            ],
            labels=[(2, "a"), (3, "b")],
        )
        newick = dinetwork_to_extended_newick(network)
        self.assertTrue(
            newick
            in [
                "(a:1.0:2.0:3.0,b:2.0:3.0:4.0);",
                "(b:2.0:3.0:4.0,a:1.0:2.0:3.0);",
            ]
        )

    def test_network(self):
        network = DiNetwork(
            edges=[(1, 2), (1, 3), (2, 4), (3, 4), (2, 5), (3, 6), (4, 7)],
            labels=[(5, "a"), (6, "b"), (7, "c")],
        )
        newick = dinetwork_to_extended_newick(network)
        self.assertTrue(
            newick
            in [
                "((a,(c)#R0),(#H0,b));",
                "(((c)#R0,a),(#H0,b));",
                "((#H0,b),(a,(c)#R0));",
                "((#H0,b),((c)#R0,a));",
                "((a,(c)#R0),(b,#H0));",
                "(((c)#R0,a),(b,#H0));",
                "((b,#H0),(a,(c)#R0));",
                "((b,#H0),((c)#R0,a));",
                "((b,(c)#R0),(#H0,a));",
                "(((c)#R0,b),(#H0,a));",
                "((#H0,a),(b,(c)#R0));",
                "((#H0,a),((c)#R0,b));",
                "((b,(c)#R0),(a,#H0));",
                "(((c)#R0,b),(a,#H0));",
                "((a,#H0),(b,(c)#R0));",
                "((a,#H0),((c)#R0,b));",
            ]
        )

    def test_network_probability(self):
        network = DiNetwork(
            edges=[(1, 2), (1, 3), (2, 4), (3, 4), (2, 5), (3, 6), (4, 7)],
            labels=[(5, "a"), (6, "b"), (7, "c")],
        )
        network[2][4]["probability"] = 0.2
        network[3][4]["probability"] = 0.8
        newick = dinetwork_to_extended_newick(network)
        self.assertTrue(
            newick
            in [
                "((b:::,(c:::)#R0:::0.8):::,(#H0:::0.2,a:::):::);",
                "(((c:::)#R0:::0.8,b:::):::,(#H0:::0.2,a:::):::);",
                "((#H0:::0.2,a:::):::,(b:::,(c:::)#R0:::0.8):::);",
                "((#H0:::0.2,a:::):::,((c:::)#R0:::0.8,b:::):::);",
                "((b:::,(c:::)#R0:::0.8):::,(a:::,#H0:::0.2):::);",
                "(((c:::)#R0:::0.8,b:::):::,(a:::,#H0:::0.2):::);",
                "((a:::,#H0:::0.2):::,(b:::,(c:::)#R0:::0.8):::);",
                "((a:::,#H0:::0.2):::,((c:::)#R0:::0.8,b:::):::);",
            ]
        )

    def test_network_all_attrs(self):
        network = DiNetwork(
            edges=[
                (
                    1,
                    2,
                    {
                        "length": 1.0,
                        "bootstrap": 2.0,
                    },
                ),
                (1, 3, {"length": 2.0, "bootstrap": 3.0, "probability": 0.2}),
                (2, 3, {"length": 3.0, "bootstrap": 4.0, "probability": 0.8}),
                (
                    2,
                    4,
                    {
                        "length": 4.0,
                        "bootstrap": 5.0,
                    },
                ),
                (
                    3,
                    5,
                    {
                        "length": 5.0,
                        "bootstrap": 6.0,
                    },
                ),
            ],
            labels=[(4, "a"), (5, "b")],
        )
        newick = dinetwork_to_extended_newick(network)
        self.assertTrue("(b:5.0:6.0:)#R0:3.0:4.0:0.8" in newick)
        self.assertTrue("a:4.0:5.0:" in newick)
        self.assertTrue("#H0:2.0:3.0:0.2" in newick)


class TestNetworkToNewickAndBack(unittest.TestCase):
    def test_larger_network(self):
        network = DiNetwork(
            edges=[
                (1, 2),
                (1, 3),
                (2, 4),
                (3, 4),
                (2, 5),
                (3, 13),
                (4, 7),
                (7, 8),
                (7, 9),
                (8, 10),
                (9, 10),
                (8, 11),
                (9, 12),
                (10, 13),
                (13, 14),
            ],
            labels=[(5, "a"), (14, "b"), (11, "c"), (12, "d")],
        )
        newick = dinetwork_to_extended_newick(network)
        network2 = extended_newick_to_dinetwork(newick)
        self.assertTrue(is_isomorphic(network, network2))
