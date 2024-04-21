import unittest

from phylox.generators.heath.heath import (
    graph_distance_to_hybridization_rate,
    suppress_degree_two_nodes,
    restrict_network_to_leaf_set,
    generate_heath_network,
)
from phylox.isomorphism import is_isomorphic


class TestHeathNetwork(unittest.TestCase):
    def test_no_time(self):
        network, hybrid_nodes, leaves, no_of_extinct = generate_heath_network(
            time_limit=0.0,
        )
        # if there is no time for any events, the network is a single root edge.
        self.assertEqual(len(network.nodes), 2)

    def test_all(self):
        network, hybrid_nodes, leaves, no_of_extinct = generate_heath_network(
            time_limit=1000000,
            taxa_limit=100,
            update_shape=2.0,
            speciation_rate_mean=4.0,
            speciation_rate_shape=2.0,
            ext_used=True,
            count_extinct=True,
            extinction_rate_mean=1.0,
            extinction_rate_shape=2.0,
            hgt_used=True,
            hgt_rate_mean=0.5,
            hgt_rate_shape=1.0,
            hgt_inheritance=0.05,
            hyb_used=True,
            hybridization_left_bound=0.05,
            hybridization_right_bound=0.1,
            hybridization_left_rate=0.1,
            hybridization_right_rate=0.05,
            simple_output=False,
        )

    # def test_seed(self):
    #     network1 = generate_network_zods(
    #         time_limit=50,
    #         speciation_rate=1,
    #         hybridization_rate=0.5,
    #         seed=1,
    #     )
    #     network2 = generate_network_zods(
    #         time_limit=50,
    #         speciation_rate=1,
    #         hybridization_rate=0.5,
    #         seed=1,
    #     )
    #     network3 = generate_network_zods(
    #         time_limit=50,
    #         speciation_rate=1,
    #         hybridization_rate=0.5,
    #         seed=2,
    #     )
    #     self.assertTrue(is_isomorphic(network1, network2))
    #     self.assertFalse(is_isomorphic(network1, network3))
