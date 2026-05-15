import unittest

from phylox import DiNetwork
from phylox.networkproperties.murepresentation import add_mu_vectors_as_attribute, add_unlabeled_mu_vectors_as_attribute
from phylox.constants import LABEL_ATTR, MUVECTOR_ATTR, MUVECTOR_UNLABELED_ATTR


class TestMuRepresentation(unittest.TestCase):
    def check_labels_in_vector(self, network):
        label_index_dict = {label: i for i, label in enumerate(sorted(network.labels.keys()))}
        for node in network.nodes:
            if LABEL_ATTR in network.nodes[node]:
                node_label = network.nodes[node][LABEL_ATTR]
                label_index = label_index_dict[node_label]
                assert network.nodes[node][MUVECTOR_ATTR][0] == node_label
                assert network.nodes[node][MUVECTOR_ATTR][label_index+1] == 1
            else:
                assert network.nodes[node][MUVECTOR_ATTR][0] == ""

    def test_tree(self):
        network = DiNetwork.from_newick("((A,B),C);")
        add_mu_vectors_as_attribute(network)
        self.check_labels_in_vector(network)

    def test_network(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")]
        )
        add_mu_vectors_as_attribute(network)
        self.check_labels_in_vector(network)
        self.assertEqual(network.nodes[1][MUVECTOR_ATTR], ("", 1, 2))


    def test_network_internal_label(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B"), (1, "I")]
        )
        add_mu_vectors_as_attribute(network)
        self.check_labels_in_vector(network)
        self.assertEqual(network.nodes[1][MUVECTOR_ATTR], ("I", 1, 2, 1))
        self.assertEqual(network.nodes[2][MUVECTOR_ATTR], ("", 1, 2, 0))


class TestMuRepresentationUnlabeled(unittest.TestCase):
    def check_labels_in_vector(self, network):
        for node in network.nodes:
            if LABEL_ATTR in network.nodes[node]:
                assert 1 in network.nodes[node][MUVECTOR_UNLABELED_ATTR]

    def test_tree(self):
        network = DiNetwork.from_newick("((A,B),C);")
        add_unlabeled_mu_vectors_as_attribute(network)
        self.check_labels_in_vector(network)

    def test_network(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B")]
        )
        add_unlabeled_mu_vectors_as_attribute(network)
        self.check_labels_in_vector(network)
        self.assertEqual(network.nodes[1][MUVECTOR_UNLABELED_ATTR], (1, 2))


    def test_network_internal_label(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "A"), (6, "B"), (1, "I")]
        )
        add_unlabeled_mu_vectors_as_attribute(network)
        self.check_labels_in_vector(network)
        self.assertEqual(network.nodes[1][MUVECTOR_UNLABELED_ATTR], (1, 2))
        self.assertEqual(network.nodes[2][MUVECTOR_UNLABELED_ATTR], (1, 2))








