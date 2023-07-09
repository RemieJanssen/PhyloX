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
    N = deepcopy(network)
    reducible_pairs = list()
    for x in N.leaves:
        reducible_pairs.extend(find_reducible_pairs_with_second(N, x))
    tree_child_sequence = list()
    while reducible_pairs:
        pair = reducible_pairs.pop()
        cherry_type = check_reducible_pair(N, *pair)
        if not cherry_type == CHERRYTYPE.NONE:
            N, _ = reduce_pair(N, *pair)
            tree_child_sequence.append(pair)
            reducible_pairs.extend(find_reducible_pairs_with_second(N, pair[1]))
            reducible_pairs.extend(find_reticulated_cherry_with_first(N, pair[1]))
    print(tree_child_sequence)
    if labels:
        tree_child_sequence = [
            (network.nodes[x][LABEL_ATTR], network.nodes[y][LABEL_ATTR])
            for x, y in tree_child_sequence
        ]
    return tree_child_sequence


def check_cherry_picking_sequence(N, cherry_picking_sequence, labels=False):
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
    if not is_tree_child(N):
        raise ValueError("N must be a tree child network")
    return check_cherry_picking_sequence(
        M, find_tree_child_sequence(N, labels=labels), labels=labels
    )
