"""
A module for generating the mu-vectors for a given DiNetwork.
The mu-vectors are added as attributes to the nodes.

author: Christopher Reichling
co-author: Remie Janssen
"""

import numpy as np
from phylox.constants import LABEL_ATTR, MUVECTOR_ATTR, MUVECTOR_UNLABELED_ATTR

def add_unlabeled_mu_vectors_as_attribute(network):
    """
    The unlabeled mu-vectors are added as a tuple(int, ..., int) consisting of num_path values.
    For each labeled node, the num_paths value belonging to their label is set to one.
    Then, starting with the leaves and going up, the number of paths to each labeled node
    is calculated for all of the nodes, by adding the mu-vector of its children to that of itself.

    The mu-vector entries are each sorted, as labels must be ignored after computing the vectors.

    The resulting vectors are stored in the node attr MUVECTOR_UNLABELED_ATTR
    and the network is modified in place.

    :param network: a DiNetwork

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.networkproperties.murepresentation import add_unlabeled_mu_vectors_as_attribute
    >>> network = DiNetwork.from_newick("((A,B),C);")
    >>> add_unlabeled_mu_vectors_as_attribute(network)
    >>> network.nodes[network.labels["A"][0]].get(MUVECTOR_UNLABELED_ATTR)
    (0, 0, 1)
    >>> network.nodes[network.labels["B"][0]].get(MUVECTOR_UNLABELED_ATTR)
    (0, 0, 1)
    >>> network.nodes[network.labels["C"][0]].get(MUVECTOR_UNLABELED_ATTR)
    (0, 0, 1)
    """
    _init_mu_representation(network)

    for node in network.nodes:
        network.nodes[node][MUVECTOR_UNLABELED_ATTR] = (int(x) for x in sorted(network.nodes[node][MUVECTOR_ATTR]))

def add_mu_vectors_as_attribute(network):
    """
    The mu-vectors are added as a tuple(str, int, ..., int) with the label of the
    node as first entry, and then the mu-vector consisting of num_path values.
    For each labeled node, the num_paths value belonging to their label is set to one.
    Then, starting with the leaves and going up, the number of paths to each labeled node
    is calculated for all of the nodes, by adding the mu-vector of its children to that of itself.

    The mu-vector entries are ordered by label.

    The resulting vectors are stored in the node attr MUVECTOR_ATTR
    and the network is modified in place.

    :param network: a DiNetwork

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.networkproperties.murepresentation import add_mu_vectors_as_attribute
    >>> network = DiNetwork.from_newick("((A,B),C);")
    >>> add_mu_vectors_as_attribute(network)
    >>> network.nodes[network.labels["A"][0]].get(MUVECTOR_ATTR)
    ('A', 1, 0, 0)
    >>> network.nodes[network.labels["B"][0]].get(MUVECTOR_ATTR)
    ('B', 0, 1, 0)
    >>> network.nodes[network.labels["C"][0]].get(MUVECTOR_ATTR)
    ('C', 0, 0, 1)
    """
    _init_mu_representation(network)

    for node in network.nodes:
        node_label = network.nodes[node].get(LABEL_ATTR, "")
        network.nodes[node][MUVECTOR_ATTR] = (node_label, *(int(x) for x in network.nodes[node][MUVECTOR_ATTR]))


def _init_mu_representation(network):
    """Sets the mu-vectors for all nodes as a np.array.
    Modifies the network in place.

    Parameters
    ----------
    network : phylox.DiNetwork
        The network to initialize the mu-represenation in

    Raises
    ------
    ValueError
        If two or more nodes have the same node label (raised in _init_mu_representation_at_labels)
    """

    _init_mu_representation_at_labels(network)

    stack = list(network.leaves)
    no_of_labels = len(network.labels)
    done = set()
    while stack:
        node = stack.pop()
        if LABEL_ATTR not in network.nodes[node]:
            network.nodes[node][MUVECTOR_ATTR] = np.zeros(no_of_labels, int)
        network.nodes[node][MUVECTOR_ATTR] += sum(network.nodes[c][MUVECTOR_ATTR] for c in network.successors(node))
        done.add(node)
        for p in network.predecessors(node):
            if all([pc in done for pc in network.successors(p)]):
                stack.append(p)


def _init_mu_representation_at_labels(network):
    """Sets the mu-vectors for all nodes with labels in `MUVECTOR_ATTR` as a np.array.
    Modifies the network in place.

    Parameters
    ----------
    network : phylox.DiNetwork
        The network to initialize the mu-represenation in

    Raises
    ------
    ValueError
        If two or more nodes have the same node label.
    """

    # labels are sorted so that we can have a tuple for the mu-vector
    label_index_dict = {label: i for i, label in enumerate(sorted(network.labels.keys()))}
    no_of_labels = len(label_index_dict)
    for label, index in label_index_dict.items():
        nodes = network.labels[label]
        if len(nodes)>1:
            raise ValueError("Cannot compute the mu-representation of a multi-labeled network.")
        if not nodes:
            # if for some reason there is a label with no nodes
            continue
        node = nodes[0]
        network.nodes[node][MUVECTOR_ATTR] = np.zeros(no_of_labels, int)
        network.nodes[node][MUVECTOR_ATTR][index] = 1