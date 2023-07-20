import random
from enum import Enum

import numpy as np

from phylox.rearrangement.move import Move, apply_move
from phylox.rearrangement.movetype import MoveType


# Pick two edges uniformly at random and add an edge between these
def random_vplu_move_at_bottom(network):
    """
    Returns a VPLU move that adds an edge between the incoming edges of two leaves.

    :param network: the network to add an edge to.
    :return: a VPLU move that adds an edge between the incoming edges of two leaves.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.movetype import MoveType
    >>> from phylox.generators.trees.well_known import generate_balanced_tree
    >>> from phylox.generators.trees.add_edges import random_vplu_move_at_bottom
    >>> tree = generate_balanced_tree(8)
    >>> move = random_vplu_move_at_bottom(tree)
    >>> move.is_type(MoveType.VPLU)
    True
    """
    leaves = list(network.leaves)
    leaf_indices = np.random.choice(range(len(leaves)), 2, replace=False)
    leaf1 = leaves[leaf_indices[0]]
    parent1 = network.parent(leaf1)
    leaf2 = leaves[leaf_indices[1]]
    parent2 = network.parent(leaf2)
    return Move(
        move_type=MoveType.VPLU,
        network=network,
        start_edge=(parent1, leaf1),
        end_edge=(parent2, leaf2),
    )


# Pick two edges uniformly at random and add an edge between these
def random_vplu_move_uniform(network):
    """
    Returns a VPLU move that adds an edge between two edges in the network.
    Two edges are chosen uniformly at random from the network.

    :param network: the network to add an edge to.
    :return: a VPLU move that adds an edge between two edges in the network.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.movetype import MoveType
    >>> from phylox.generators.trees.well_known import generate_balanced_tree
    >>> from phylox.generators.trees.add_edges import random_vplu_move_uniform
    >>> tree = generate_balanced_tree(8)
    >>> move = random_vplu_move_uniform(tree)
    >>> move.is_type(MoveType.VPLU)
    True
    """
    edges = list(network.edges())
    edge_indices = np.random.choice(range(len(edges)), 2, replace=False)
    return Move(
        move_type=MoveType.VPLU,
        network=network,
        start_edge=edges[edge_indices[0]],
        end_edge=edges[edge_indices[1]],
    )


def random_vplu_move_local(network, stop_prob=0.2, max_steps=None, max_tries=None):
    """
    Returns a VPLU move that adds an edge between two edges in the network.
    Pick one edge, move a random number of edges through the network to find a second edge.

    :param network: the network to add an edge to.
    :param stop_prob: the probability to stop the random walk.
    :param max_steps: the maximum number of steps to take in the random walk.
    :param max_tries: the maximum number of tries to find a second edge.
    :return: a VPLU move that adds an edge between two edges in the network.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.movetype import MoveType
    >>> from phylox.generators.trees.well_known import generate_balanced_tree
    >>> from phylox.generators.trees.add_edges import random_vplu_move_local
    >>> tree = generate_balanced_tree(8)
    >>> move = random_vplu_move_local(tree)
    >>> move.is_type(MoveType.VPLU)
    True
    """
    try_number = 1
    while max_tries == None or try_number <= max_tries:
        # Pick a random edge
        edge1 = random.choice(list(network.edges()))
        edge2 = None
        # Initiate the random walk, by choosing an orientation
        previous_node = random.choice(edge1)
        current_node = edge1[0]
        if current_node == previous_node:
            current_node = edge1[1]
        # Take a number of steps
        step_number = 1
        while max_steps == None or step_number <= max_steps:
            previous_node, current_node = current_node, random.choice(
                list(network.successors(current_node))
                + list(network.predecessors(current_node))
            )
            if random.random() < stop_prob:
                break
            step_number += 1
        # Set the new edge
        edge2 = (previous_node, current_node)
        if edge2 not in network.edges():
            edge2 = (current_node, previous_node)
        # Add an edge if possible, otherwise repeat the search
        if edge1 != edge2:
            break
        try_number += 1
    return Move(
        move_type=MoveType.VPLU,
        network=network,
        start_edge=edge1,
        end_edge=edge2,
    )


class AddEdgeMethod(Enum):
    UNIFORM = 1
    BOTTOM = 2
    LOCAL = 3


def network_from_tree(tree, reticulations, method):
    """
    Returns a network with the given number of reticulations added to the given tree.

    :param tree: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute LABEL_ATTR.
    :param reticulations: the number of reticulations to add to the tree.
    :param method: the method to use to add the reticulations.
    :return: a network with the given number of reticulations added to the given tree.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.movetype import MoveType
    >>> from phylox.generators.trees.well_known import generate_balanced_tree
    >>> from phylox.generators.trees.add_edges import network_from_tree, AddEdgeMethod
    >>> tree = generate_balanced_tree(8)
    >>> network = network_from_tree(tree, 2, AddEdgeMethod.BOTTOM)
    >>> len(network.leaves)
    8
    >>> network.reticulation_number
    2
    """

    if method == AddEdgeMethod.BOTTOM:
        leaves = tree.leaves

    if method == AddEdgeMethod.UNIFORM:
        add_edge_method = random_vplu_move_uniform
    elif method == AddEdgeMethod.BOTTOM:
        add_edge_method = random_vplu_move_at_bottom
    elif method == AddEdgeMethod.LOCAL:
        add_edge_method = random_vplu_move_local
    else:
        raise ValueError(f"Method {method} not implemented, use one of {AddEdgeMethod}")

    network = tree.copy()
    reticulations = int(reticulations)
    while reticulations > 0:
        try:
            move = add_edge_method(network)
            network = apply_move(network, move)
            reticulations -= 1
        except:
            pass
    return network
