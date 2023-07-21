import math
import random

import networkx as nx
import numpy as np

from phylox import DiNetwork


def generate_balanced_tree(number_of_leaves):
    """
    Generates a balanced tree with the given number of leaves.

    :param number_of_leaves: the number of leaves of the balanced tree (power of 2).
    :return: a balanced tree with the given number of leaves.

    :example:
    >>> from phylox.generators.trees.well_known import generate_balanced_tree
    >>> tree = generate_balanced_tree(8)
    >>> len(tree.leaves)
    8
    """
    tree = DiNetwork()
    leaves = list(range(1, number_of_leaves + 1))
    tree.add_nodes_from(leaves)
    roots = leaves
    current_node = number_of_leaves + 1
    while len(roots) > 1:
        new_roots = []
        while len(roots) > 1:
            root1 = roots.pop()
            root2 = roots.pop()
            tree.add_edges_from([(current_node, root1), (current_node, root2)])
            new_roots.append(current_node)
            current_node += 1
        roots = new_roots
    tree.add_edge(0, roots[0])
    return tree


def generate_caterpillar(number_of_leaves):
    """
    Generates a caterpillar tree with the given number of leaves.

    :param number_of_leaves: the number of leaves of the caterpillar tree.
    :return: a caterpillar tree with the given number of leaves.

    :example:
    >>> from phylox.generators.trees.well_known import generate_caterpillar
    >>> tree = generate_caterpillar(5)
    >>> len(tree.leaves)
    5
    """
    tree = DiNetwork()
    tree.add_node(number_of_leaves + 1)
    for i in range(1, number_of_leaves):
        tree.add_edges_from(
            [
                (number_of_leaves + i, i),
                (number_of_leaves + i, number_of_leaves + i + 1),
            ]
        )
    tree = nx.relabel_nodes(tree, {2 * number_of_leaves: number_of_leaves})
    tree.add_edge(0, number_of_leaves + 1)
    return tree
