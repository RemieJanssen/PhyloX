"""
A module for checking isomorphism between phylogenetic networks and counting automorphisms of phylogenetic networks.

Currently mostly uses the networkx isomorphism checker, extended with some specifics for phylogenetic networks 
such as labels and partial isomorphisms.
"""

from copy import deepcopy

import networkx as nx

from phylox.constants import LABEL_ATTR

#: The node attribute used to store the isometry label of a node.
ISOMETRY_LABEL_ATTR = "isometry_label"
#: The prefix used for isometry labels.
ISOMETRY_LABEL_PREFIX = "isometry_label_prefix_"
#: The prefix used for automorphism labels.
AUTOMORPHISM_LABEL_PREFIX = "automorphism_label_prefix_"

# Checks whether the nodes with the given attributes have the same label
def _same_isometry_labels(node1_attributes, node2_attributes):
    """
    Checks whether two nodes have the same label

    :param node1_attributes: the attributes of a node
    :param node2_attributes: the attributes of a node
    :return: True if the isometry label attribute ISOMETRY_LABEL_ATTR is the same, False otherwise.
    """
    return node1_attributes.get(ISOMETRY_LABEL_ATTR) == node2_attributes.get(
        ISOMETRY_LABEL_ATTR
    )


# Checks whether the nodes with the given attributes have the same label
def _same_isometry_labels_and_labels(node1_attributes, node2_attributes):
    """
    Checks whether two nodes have the same label

    :param node1_attributes: the attributes of a node
    :param node2_attributes: the attributes of a node
    :return: True if the isometry label attribute ISOMETRY_LABEL_ATTR is the same, False otherwise.
    """
    return node1_attributes.get(ISOMETRY_LABEL_ATTR) == node2_attributes.get(
        ISOMETRY_LABEL_ATTR
    ) and node1_attributes.get(LABEL_ATTR) == node2_attributes.get(LABEL_ATTR)


# Checks whether two networks are labeled isomorpgic
def is_isomorphic(network1, network2, partial_isomorphism=None, ignore_labels=False):
    """
    Determines whether two networks are labeled isomorphic.

    :param network1: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute LABEL_ATTR.
    :param network2: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute LABEL_ATTR.
    :return: True if the networks are labeled isomorphic, False otherwise.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.isomorphism.base import is_isomorphic
    >>> network1 = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ...     labels=[(4, "A"), (5, "B")],
    ... )
    >>> network2 = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,5),(3,6)],
    ...     labels=[(5, "B"), (6, "A")],
    ... )
    >>> is_isomorphic(network1, network2, ignore_labels=True)
    True
    >>> is_isomorphic(network1, network2, ignore_labels=False)
    False
    >>> is_isomorphic(network1, network2, partial_isomorphism=[(4,6)], ignore_labels=True)
    False
    """
    nw1 = deepcopy(network1)
    nw2 = deepcopy(network2)

    same_labels = _same_isometry_labels_and_labels
    if ignore_labels:
        same_labels = _same_isometry_labels

    partial_isomorphism = partial_isomorphism or []
    for i, corr in enumerate(partial_isomorphism):
        if not same_labels(nw1.nodes[corr[0]], nw2.nodes[corr[1]]):
            return False
        nw1.nodes[corr[0]][ISOMETRY_LABEL_ATTR] = f"{ISOMETRY_LABEL_PREFIX}{i}"
        nw2.nodes[corr[1]][ISOMETRY_LABEL_ATTR] = f"{ISOMETRY_LABEL_PREFIX}{i}"

    return nx.is_isomorphic(nw1, nw2, node_match=same_labels)


def _count_automorphisms(
    network,
    ignore_labels=False,
    partial_isomorphism=None,
    nodes_available=None,
    nodes_to_do=None,
):
    """
    Determines the number of automorphisms of a network.

    An implementation of the algorithm described in
    "A note on the graph isomorphism counting problem" by
    Rudolf Matheson in 1979.

    :param network: a phylox.DiNetwork phylogenetic network, i.e., a DAG with leaf labels.
    :param ignore_labels: if True, the automorphisms are counted without considering the labels of the nodes.
    :param partial_isomorphism: a partial isomorphism between the network and itself.
    :param nodes_available: the nodes that are available to be matched.
    :param nodes_to_do: the nodes that still need to be matched.
    :return: the number of automorphisms of the network.
    """
    nodes_available = nodes_available or []
    nodes_to_do = nodes_to_do if nodes_to_do is not None else set(network.nodes())
    same_labels = _same_isometry_labels_and_labels
    if ignore_labels:
        same_labels = _same_isometry_labels

    number_of_automorphisms = 1
    while nodes_to_do:
        node_pair_to_remove = partial_isomorphism.pop()
        node_to_remove = node_pair_to_remove[0]
        nodes_to_do.remove(node_to_remove)
        matches = 1
        for try_to_match_node in nodes_available:
            if is_isomorphic(
                network,
                network,
                partial_isomorphism=partial_isomorphism
                + [(node_to_remove, try_to_match_node)],
                ignore_labels=ignore_labels,
            ):
                matches += 1
        number_of_automorphisms *= matches
        nodes_available.append(node_to_remove)
    return number_of_automorphisms


def count_automorphisms(network, ignore_labels=False):
    """
    Determines the number of automorphisms of a network.

    :param network: a phylogenetic network, i.e., a DAG with leaf labels.
    :param ignore_labels: if True, the automorphisms are counted without considering the labels of the nodes.
    :return: the number of automorphisms of the network.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.isomorphism.base import count_automorphisms
    >>> network = DiNetwork(
    ...     edges=[(-1,0), (0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)],
    ...     labels=[(3, "A"), (4, "B"), (5, "C"), (6, "C")],
    ... )
    >>> count_automorphisms(network, ignore_labels=True)
    8
    >>> count_automorphisms(network, ignore_labels=False)
    2
    """
    partial_isomorphism = [(a, a) for a in network.nodes()]
    return _count_automorphisms(
        network,
        ignore_labels=ignore_labels,
        partial_isomorphism=partial_isomorphism,
    )
