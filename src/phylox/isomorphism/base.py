import networkx as nx
from phylox.dinetwork import LABEL_ATTR
from copy import deepcopy

# Checks whether the nodes with the given attributes have the same label
def same_labels(node1_attributes, node2_attributes):
    """
    Checks whether two nodes have the same label

    :param node1_attributes: the attributes of a node
    :param node2_attributes: the attributes of a node
    :return: True if the label attribute is the same, False otherwise.
    """
    return node1_attributes.get('label') == node2_attributes.get('label')


# Checks whether two networks are labeled isomorpgic
def is_isomorphic(network1, network2, partial_isomorphism=None):
    """
    Determines whether two networks are labeled isomorphic.

    :param network1: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute `label'.
    :param network2: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute `label'.
    :return: True if the networks are labeled isomorphic, False otherwise.
    """
    nw1 = deepcopy(network1)
    nw2 = deepcopy(network2)

    partial_isomorphism = partial_isomorphism or []
    for i, corr in enumerate(partial_isomorphism):
        if not same_labels(nw1.nodes[corr[0]], nw2.nodes[corr[1]]):
            return False
        nw1.nodes[corr[0]][LABEL_ATTR] = f"{i}_isom_label"
        nw2.nodes[corr[1]][LABEL_ATTR] = f"{i}_isom_label"

    return nx.is_isomorphic(nw1, nw2, node_match=same_labels)
