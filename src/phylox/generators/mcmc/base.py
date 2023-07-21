import os
import sys

import networkx as nx
import numpy as np

from phylox.base import find_unused_node
from phylox.exceptions import InvalidMoveDefinitionException, InvalidMoveException
from phylox.isomorphism import count_automorphisms
from phylox.rearrangement.invertsequence import from_edge
from phylox.rearrangement.move import Move, apply_move
from phylox.rearrangement.movetype import MoveType


def acceptance_probability(
    network,
    result_network,
    move,
    move_type_probabilities,
    number_of_leaves=None,
    current_reticulation_number=None,
    symmetries=False,
):
    """
    Computes the acceptance probability of a move.

    :param network: the network before the move.
    :param result_network: the network after the move.
    :param move: the move.
    :param move_type_probabilities: the move type probabilities.
    :param number_of_leaves: the number of leaves in the network.
    :param current_reticulation_number: the current number of reticulations in the network.
    :param symmetries: whether to correct for symmetries.
    :return: the acceptance probability of the move.
    """
    current_reticulation_number = (
        current_reticulation_number or network.reticulation_number
    )
    number_of_leaves = number_of_leaves or len(network.leaves)
    p = 0
    if move.move_type in [MoveType.TAIL, MoveType.HEAD]:
        p = 1
    if move.move_type == MoveType.VPLU:
        no_edges_network = float(
            2 * number_of_leaves + 3 * current_reticulation_number - 1
        )
        no_edges_network_after = no_edges_network + 3
        p = (
            (
                move_type_probabilities[MoveType.VMIN]
                / move_type_probabilities[MoveType.VPLU]
            )
            * no_edges_network**2
            / (no_edges_network_after)
        )
    if move.move_type == MoveType.VMIN:
        no_edges_network = float(
            2 * number_of_leaves + 3 * current_reticulation_number - 1
        )
        no_edges_network_after = no_edges_network - 3
        if no_edges_network > 3:
            p = (
                (
                    move_type_probabilities[MoveType.VPLU]
                    / move_type_probabilities[MoveType.VMIN]
                )
                * no_edges_network
                / (no_edges_network_after**2)
            )
    if symmetries:
        # correct for number of representations, i.e., symmetries.
        p *= count_automorphisms(network) / count_automorphisms(result_network)
    return p


def sample_mcmc_networks(
    starting_network,
    move_type_probabilities,
    restriction_map=None,
    correct_symmetries=True,
    burn_in=1000,
    number_of_samples=1,
    add_root_if_necessary=False,
):
    """
    Samples phylogenetic networks using a Markov-Chain Monte Carlo method.

    :param starting_network: the phylox.DiNetwork used as the starting point of the Markov chain.
    :param move_type_probabilities: a dictionary mapping MoveTypes to probabilities.
    :param restriction_map: a boolean function that takes a phylox.DiNetwork as input.
    :param correct_symmetries: whether to correct for symmetries in the acceptance probability, set to True for uniform distribution.
    :param burn_in: the number of steps (including rejected proposals) between each sample.
    :param number_of_samples: the number of networks to sample.
    :param add_root_if_necessary: whether to add a root edge to each root if it has out-degree > 1.

    :return: a list of phylox.DiNetwork objects.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.movetype import MoveType
    >>> from phylox.generators.mcmc import sample_mcmc_networks
    >>> starting_network = DiNetwork(
    ...     edges = ((0,1), (0,2), (1,2), (1,3), (2,4)),
    ...     labels = ((3, "A"), (4, "B")),
    ... )
    >>> move_type_probabilities = {
    ...     MoveType.TAIL: 0.2,
    ...     MoveType.VPLU: 0.4,
    ...     MoveType.VMIN: 0.4,
    ... }
    >>> restriction_map = (lambda nw: nw.reticulation_number < 2)
    >>> sampled_networks = sample_mcmc_networks(
    ...     starting_network,
    ...     move_type_probabilities,
    ...     restriction_map=restriction_map,
    ...     correct_symmetries=False,
    ...     burn_in=100,
    ...     number_of_samples=50,
    ...     add_root_if_necessary=False,
    ... )
    >>> all([network.reticulation_number<2 for network in sampled_networks])
    True
    >>> all([len(network.leaves)==2 for network in sampled_networks])
    True
    """
    network = starting_network.copy()
    current_reticulation_number = network.reticulation_number
    number_of_leaves = len(network.leaves)
    if add_root_if_necessary:
        for root in network.roots:
            if network.out_degree(root) >= 1:
                new_root = find_unused_node(network)
                network.add_edges_from([(new_root, root)])
                root = new_root
        roots = network._set_roots()
    available_reticulations = set()
    available_tree_nodes = set()

    sample = []

    for _ in range(number_of_samples):
        non_moves = 0
        for _ in range(burn_in):
            try:
                move = Move.random_move(
                    network,
                    available_tree_nodes=available_tree_nodes,
                    available_reticulations=available_reticulations,
                    move_type_probabilities=move_type_probabilities,
                )
                result_network = apply_move(network, move)
            except (InvalidMoveException, InvalidMoveDefinitionException) as e:
                non_moves += 1
                continue
            if np.random.random() > acceptance_probability(
                network,
                result_network,
                move,
                move_type_probabilities,
                number_of_leaves=number_of_leaves,
                current_reticulation_number=current_reticulation_number,
                symmetries=correct_symmetries,
            ):
                non_moves += 1
                continue
            # only apply the move if the restrinction_map returns True
            if not (restriction_map is None or restriction_map(result_network)):
                non_moves += 1
                continue
            if move.move_type in [MoveType.TAIL, MoveType.HEAD]:
                network = result_network
            if move.move_type == MoveType.VPLU:
                current_reticulation_number += 1
                available_tree_nodes.discard(move.start_node)
                available_reticulations.discard(move.end_node)
                network = result_network
            if move.move_type == MoveType.VMIN:
                current_reticulation_number -= 1
                available_tree_nodes.add(move.removed_edge[0])
                available_reticulations.add(move.removed_edge[1])
                network = result_network
        sample.append(network)
    return sample
