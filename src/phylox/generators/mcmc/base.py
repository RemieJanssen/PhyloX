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

    for index in range(number_of_samples):
        non_moves = 0
        for j in range(burn_in):
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
