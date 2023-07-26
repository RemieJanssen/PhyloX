import unittest

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR
from phylox.parser import extended_newick_to_dinetwork
from phylox.isomorphism import is_isomorphic


class TestExtendedNewickToDiNetwork(unittest.TestCase):
    def test_small_tree(self):
        newick = "(a,b,(c,d));"
        network = extended_newick_to_dinetwork(newick)
        network2 = DiNetwork(
            edges=[(1, 2), (1, 3), (1, 4), (4, 5), (4, 6)],
            labels=[(2, "a"), (3, "b"), (5, "c"), (6, "d")],
        )
        print(network.edges)
        print(network2.edges)
        print(network.nodes(data=True))
        print(network2.nodes(data=True))

        self.assertTrue(is_isomorphic(network, network2))
