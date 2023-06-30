import random
from enum import Enum

import numpy as np

from phylox.rearrangement.move import Move, apply_move
from phylox.rearrangement.movetype import MoveType


# Pick two edges uniformly at random and add an edge between these
def random_vplu_move_at_bottom(network):
    """
    Returns a VPLU move that adds an edge between the incoming edges of two leaves.
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
