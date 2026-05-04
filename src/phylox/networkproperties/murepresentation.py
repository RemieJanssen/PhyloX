"""
A module for generating the mu-vectors for a given DiNetwork.
The mu-vectors are added as attributes to the nodes.
author: Christopher Reichling
"""

import networkx as nx
import numpy as np
from phylox.constants import LABEL_ATTR, MUVECTOR_ATTR

def add_mu_vectors_as_attribute(network):

    """
    The mu-vectors are added as a dictionary containing leaf_label:num_paths pairs.
    For each leaf the num_paths value belonging to their leaf_label is set to one.
    Then, starting with the leaves, for each parent of each node the num_paths value
    corresponding to the leaf_label is increased by 1 recursively.

    :param network: a DiNetwork
    """

    for node in network.nodes:
        # network.nodes[node][MUVECTOR_ATTR] = np.zeros(len(network.leaves))
        network.nodes[node][MUVECTOR_ATTR] = {}
        for leaf in network.leaves:
            label = network.nodes[leaf][LABEL_ATTR]
            network.nodes[node][MUVECTOR_ATTR][label] = 0

    for leaf in network.leaves:
        label = network.nodes[leaf][LABEL_ATTR]
        network.nodes[leaf][MUVECTOR_ATTR][label] = 1
        _update_parents(network, leaf, label)


def _update_parents(network, node, label):
    for parent in network.predecessors(node):
        network.nodes[parent][MUVECTOR_ATTR][label] += 1
        _update_parents(network,parent,label)

