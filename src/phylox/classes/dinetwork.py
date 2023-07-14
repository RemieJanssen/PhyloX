from copy import deepcopy

import networkx as nx

from phylox.cherrypicking import CHERRYTYPE, is_second_in_reducible_pair, reduce_pair
from phylox.constants import LABEL_ATTR


def is_binary(network):
    """
    Checks if the network is binary.

    :param network: a phylogenetic network phylox.DiNetwork.
    :return: true if the network is binary, false otherwise.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.classes.dinetwork import is_binary
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ... )
    >>> is_binary(network)
    True

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(1,4)],
    ... )
    >>> is_binary(network)
    False
    """

    binary_node_types = [
        [0, 1],  # root
        [0, 2],  # root
        [1, 2],  # tree node
        [2, 1],  # reticulation
        [1, 0],  # leaf
    ]
    for node in network.nodes:
        degrees = [network.in_degree(node), network.out_degree(node)]
        if degrees not in binary_node_types:
            return False
    return True


def is_orchard(network):
    """
    Checks if the network is an orchard.

    :param network: a phylogenetic network phylox.DiNetwork.
    :return: true if the network is an orchard, false otherwise.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.classes.dinetwork import is_orchard
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ... )
    >>> is_orchard(network)
    True

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,4),(3,5),(2,5),(3,4),(4,6),(5,7)],
    ... )
    >>> is_orchard(network)
    False
    """
    if len(network) == 0:
        return True
    leaves = network.leaves
    root = list(network.roots)[0]

    # make a copy and fix a root edge
    network_copy = deepcopy(network)
    if network_copy.out_degree(root) > 1:
        new_node = -1
        while new_node in network_copy.nodes:
            new_node -= 1
        network_copy.add_edge(new_node, root)

    # try to reduce the network copy
    done = False
    while not done:
        checked_all_leaves = True
        for leaf in leaves:
            pair = is_second_in_reducible_pair(network_copy, leaf)
            if pair:
                network_copy, cherry_type = reduce_pair(network_copy, *pair)
                if cherry_type == CHERRYTYPE.CHERRY:
                    leaves.remove(pair[0])
                checked_all_leaves = False
                break
        if len(network_copy.edges) == 1:
            return True
        done = checked_all_leaves
    return False


def is_stack_free(network):
    """
    Checks if the network is stack-free.

    :param network: a phylogenetic network phylox.DiNetwork.
    :return: true if the network is stack-free, false otherwise.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.classes.dinetwork import is_stack_free
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,4),(3,5),(2,5),(3,4),(4,6),(5,7)],
    ... )
    >>> is_stack_free(network)
    True

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(3,5),(2,4),(4,5),(4,6),(5,7)],
    ... )
    >>> is_stack_free(network)
    False
    """
    for node in network.nodes:
        if network.is_reticulation(node) and any(
            [network.is_reticulation(child) for child in network.successors(node)]
        ):
            return False
    return True


def _is_endpoint_of_w_fence(network, node):
    if not network.is_reticulation(node):
        return False
    previous_node = node
    current_node = network.child(node)
    currently_at_fence_top = False
    while True:
        if network.is_leaf(current_node):
            return False
        if network.is_reticulation(current_node):
            if currently_at_fence_top:
                return True
            next_node = network.parent(current_node, exclude=[previous_node])
        if network.is_tree_node(current_node):
            if not currently_at_fence_top:
                return False
            next_node = network.child(current_node, exclude=[previous_node])
        previous_node, current_node = current_node, next_node
        currently_at_fence_top = not currently_at_fence_top


def is_tree_based(network):
    """
    Checks if the network is tree-based.

    :param network: a phylogenetic network phylox.DiNetwork.
    :return: true if the network is tree-based, false otherwise.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.classes.dinetwork import is_tree_based
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ... )
    >>> is_tree_based(network)
    True
    """
    if not is_binary(network):
        raise NotImplementedError(
            "tree-basedness cannot be computed for non-binary networks yet."
        )

    if len(network) > 0 and not nx.is_weakly_connected(network):
        return False

    if len(network.roots) > 1:
        return False

    for node in network.nodes:
        if _is_endpoint_of_w_fence(network, node):
            return False
    return True


def is_tree_child(network):
    """
    Checks if the network is a tree-child network.

    :param network: a phylogenetic network phylox.DiNetwork.
    :return: true if the network is a tree-child network, false otherwise.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.classes.dinetwork import is_tree_child
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ... )
    >>> is_tree_child(network)
    True

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,4),(3,5),(2,5),(3,4),(4,6),(5,7)],
    ... )
    >>> is_tree_child(network)
    False

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(3,5),(2,4),(4,5),(4,6),(5,7)],
    ... )
    >>> is_tree_child(network)
    False
    """

    for node in network.nodes:
        if network.is_leaf(node):
            continue
        if all([network.is_reticulation(child) for child in network.successors(node)]):
            return False
    return True


def is_leaf_labeled_single_root_network(network):
    """
    Checks if the network is a leaf-labeled network with a single root.

    :param network: a phylogenetic network phylox.DiNetwork.
    :return: a boolean value.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.classes.dinetwork import is_leaf_labeled_single_root_network
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ...     labels=[(4, "A"), (5, "B")],
    ... )
    >>> is_leaf_labeled_single_root_network(network)
    True

    >>> network = DiNetwork(
    ...     edges=[(0,1),(2,3)],
    ...     labels=[(1, "A"), (3, "B")],
    ... )
    >>> is_leaf_labeled_single_root_network(network)
    False

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3)],
    ...     labels=[(3, "B")],
    ... )
    >>> is_leaf_labeled_single_root_network(network)
    False
    """
    rootFound = False
    for v in network.nodes:
        if network.in_degree(v) == 0:
            if rootFound:
                return False
            rootFound = True
        if network.out_degree(v) == 0 and not LABEL_ATTR in network.nodes[v]:
            return False
        if network.out_degree(v) == 1 and network.in_degree(v) == 1:
            return False
        if network.out_degree(v) > 1 and network.in_degree(v) > 1:
            return False
    return True
