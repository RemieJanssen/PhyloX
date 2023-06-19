from copy import deepcopy

import networkx as nx

from phylox.dinetwork import LABEL_ATTR

ISOMETRY_LABEL_ATTR = "isometry_label"
ISOMETRY_LABEL_TAG = "isometry_label_tag_"
AUTOMORPHISM_LABEL_TAG = "automorphism_label_tag_"

# Checks whether the nodes with the given attributes have the same label
def same_isometry_labels(node1_attributes, node2_attributes):
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
def same_isometry_labels_and_labels(node1_attributes, node2_attributes):
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

    :param network1: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute `label'.
    :param network2: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute `label'.
    :return: True if the networks are labeled isomorphic, False otherwise.
    """
    nw1 = deepcopy(network1)
    nw2 = deepcopy(network2)

    same_labels = same_isometry_labels_and_labels
    if ignore_labels:
        same_labels = same_isometry_labels

    partial_isomorphism = partial_isomorphism or []
    for i, corr in enumerate(partial_isomorphism):
        if not same_labels(nw1.nodes[corr[0]], nw2.nodes[corr[1]]):
            return False
        nw1.nodes[corr[0]][ISOMETRY_LABEL_ATTR] = f"{ISOMETRY_LABEL_TAG}{i}"
        nw2.nodes[corr[1]][ISOMETRY_LABEL_ATTR] = f"{ISOMETRY_LABEL_TAG}{i}"

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

    :param network: a phylogenetic network, i.e., a DAG with leaf labels.
    :param ignore_labels: if True, the automorphisms are counted without considering the labels of the nodes.
    :param partial_isomorphism: a partial isomorphism between the network and itself.
    :param nodes_available: the nodes that are available to be matched.
    :param nodes_to_do: the nodes that still need to be matched.
    :return: the number of automorphisms of the network.

    """
    nodes_available = nodes_available or []
    nodes_to_do = nodes_to_do if nodes_to_do is not None else set(network.nodes())
    same_labels = same_isometry_labels_and_labels
    if ignore_labels:
        same_labels = same_isometry_labels

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
    """
    partial_isomorphism = [(a, a) for a in network.nodes()]
    return _count_automorphisms(
        network,
        ignore_labels=ignore_labels,
        partial_isomorphism=partial_isomorphism,
    )
