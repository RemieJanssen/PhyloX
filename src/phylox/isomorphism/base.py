import networkx as nx
from phylox.dinetwork import LABEL_ATTR
from copy import deepcopy

ISOMETRY_LABEL_ATTR = 'isometry_label'
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
    return node1_attributes.get(ISOMETRY_LABEL_ATTR) == node2_attributes.get(ISOMETRY_LABEL_ATTR)


# Checks whether the nodes with the given attributes have the same label
def same_isometry_labels_and_labels(node1_attributes, node2_attributes):
    """
    Checks whether two nodes have the same label

    :param node1_attributes: the attributes of a node
    :param node2_attributes: the attributes of a node
    :return: True if the isometry label attribute ISOMETRY_LABEL_ATTR is the same, False otherwise.
    """
    return node1_attributes.get(ISOMETRY_LABEL_ATTR) == node2_attributes.get(ISOMETRY_LABEL_ATTR) and \
           node1_attributes.get(LABEL_ATTR) == node2_attributes.get(LABEL_ATTR)


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



def _count_automorphisms(network, ignore_labels=False, partial_isomorphism=None):

    pass

def count_automorphisms(network, ignore_labels=False):
    pass