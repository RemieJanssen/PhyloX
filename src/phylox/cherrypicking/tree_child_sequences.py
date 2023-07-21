"""
Functions for finding tree child sequences in a network.

Based on "Linear Time Algorithm for Tree-Child network Containment" by Remie Janssen and Yukihiro Murakami (2020)
"""

import ast
import itertools
import random
import time
from copy import deepcopy

import networkx as nx

from phylox import DiNetwork
from phylox.cherrypicking import (
    CHERRYTYPE,
    check_reducible_pair,
    find_reducible_pairs_with_second,
    find_reticulated_cherry_with_first,
    reduce_pair,
)
from phylox.classes.dinetwork import is_tree_child
from phylox.constants import LABEL_ATTR


def find_tree_child_sequence(network, labels=False):
    """
    Find a tree child sequence for a network. If the network is orchard but not tree-child,
    this function may produce a cherry-picking sequence that reduces the network to a single edge,
    but that is not a tree-child sequence.

    :param network: The network to find a tree child sequence for.
    :type network: phylox.DiNetwork
    :param labels: Whether to return the tree child sequence as a list of labels instead of a list of nodes.
    :type labels: bool
    :return: The tree child sequence.
    :rtype: list

    :example:
    >>> from phylox import DiNetwork
    >>> N = DiNetwork(
    ...     edges=[(-1, 0), (0, 1), (0, 2), (1, 3), (1, 4), (2, 3), (2, 4), (3, 5), (4, 6)],
    ...     labels=[(5, "A"), (6, "B")],
    ... )
    >>> try:
    ...     find_tree_child_sequence(N, labels=True)
    ... except ValueError:
    ...     exception_raised = True
    >>> exception_raised
    True

    >>> N = DiNetwork(
    ...     edges=[(-1, 0), (0, 1), (0, 2), (1, 2), (1, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
    ...     labels=[(5, "A"), (6, "B")],
    ... )
    >>> TC_sequence = find_tree_child_sequence(N, labels=True)
    >>> TC_sequence == [("B", "A"), ("B", "A"), ("B", "A")] or TC_sequence == [("A", "B"), ("B", "A"), ("B", "A")]
    True
    """
    N = deepcopy(network)
    reducible_pairs = []
    for x in N.leaves:
        reducible_pairs.extend(find_reducible_pairs_with_second(N, x))
    tree_child_sequence = list()
    while reducible_pairs:
        pair = reducible_pairs.pop()
        cherry_type = check_reducible_pair(N, *pair)
        N, cherry_type = reduce_pair(N, *pair)
        if cherry_type != CHERRYTYPE.NONE:
            tree_child_sequence.append(pair)
            reducible_pairs.extend(find_reducible_pairs_with_second(N, pair[1]))
            reducible_pairs.extend(find_reticulated_cherry_with_first(N, pair[1]))
    if labels:
        tree_child_sequence = [
            (network.nodes[x][LABEL_ATTR], network.nodes[y][LABEL_ATTR])
            for x, y in tree_child_sequence
        ]
    if N.size() == 1:
        return tree_child_sequence
    raise ValueError("N could not be reduced to a single edge")


def check_cherry_picking_sequence(N, cherry_picking_sequence, labels=False):
    """
    Check if a cherry picking sequence is valid for a network.

    :param N: The network to check the cherry picking sequence for.
    :type N: phylox.DiNetwork
    :param cherry_picking_sequence: The cherry picking sequence to check.
    :type cherry_picking_sequence: list
    :param labels: Whether the cherry picking sequence is given as a list of labels instead of a list of nodes.
    :type labels: bool
    :return: True if the cherry picking sequence reduces the network to a single edge, False otherwise.
    :rtype: bool

    :example:
    >>> from phylox import DiNetwork
    >>> N = DiNetwork(
    ...     edges=[(-1, 0), (0, 1), (0, 2), (1, 2), (1, 3), (2, 4)],
    ...     labels=[(3, "A"), (4, "B")],
    ... )
    >>> check_cherry_picking_sequence(N, [("B", "A"), ("B", "A")], labels=True)
    True
    """
    if labels:
        cherry_picking_sequence = [
            (N.label_to_node_dict.get(x, None), N.label_to_node_dict.get(y, None))
            for x, y in cherry_picking_sequence
        ]
    for pair in cherry_picking_sequence:
        N, _ = reduce_pair(N, *pair)
    if N.size() == 1:
        return True
    return False


def tree_child_network_contains(N, M, labels=False):
    """
    Check if a tree child network N contains another network M.

    :param N: The network to check for containment.
    :type N: phylox.DiNetwork
    :param M: The network to check for.
    :type M: phylox.DiNetwork
    :param labels: Whether to return the tree child sequence as a list of labels instead of a list of nodes.
    :type labels: bool
    :return: True if N contains M, False otherwise.
    :rtype: bool

    :example:
    >>> from phylox import DiNetwork
    >>> N = DiNetwork(
    ...     edges=[(-1, 0), (0, 1), (0, 2), (1, 2), (1, 3), (2, 4)],
    ...     labels=[(3, "A"), (4, "B")],
    ... )
    >>> M = DiNetwork(
    ...     edges=[(-1, 0), (0, 1), (0, 2)],
    ...     labels=[(1, "A"), (2, "B")],
    ... )
    >>> tree_child_network_contains(N, M, labels=True)
    True
    """
    if not is_tree_child(N):
        raise ValueError("N must be a tree child network")
    return check_cherry_picking_sequence(
        M, find_tree_child_sequence(N, labels=labels), labels=labels
    )
