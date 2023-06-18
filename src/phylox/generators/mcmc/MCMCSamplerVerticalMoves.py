import networkx as nx
from RearrDistance_Tools import *
import os
import sys
from phylox.isomorphism import count_automorphisms
from phylox.rearrangement.movetype import MoveType
from phylox.rearrangement.move import apply_move, Move
from phylox.rearrangement.invertsequence import from_edge


def acceptance_probability(
    network,
    move,
    move_type_probabilities,
    number_of_leaves=None,
    current_reticulation_number=None,
    symmetries=False,
    result_network=None,
):
    current_reticulation_number = current_reticulation_number or network.reticulation_number()
    number_of_leaves = number_of_leaves or network.leaf_number()
    p = 0
    if move.move_type in [MoveType.TAIL, MoveType.HEAD]:
        p = 1
    if move.move_type == MoveType.VPLU:
        no_edges_network = float(2 * number_of_leaves + 3 * current_reticulation_number - 1)
        no_edges_network_after = no_edges_network + 3
        p = (
            (move_type_probabilities[MoveType.VMIN] / move_type_probabilities[MoveType.VPLU])
            * no_edges_network**2
            / (no_edges_network_after)
        )
    if move.move_type == MoveType.VMIN:
        no_edges_network = float(2 * number_of_leaves + 3 * current_reticulation_number - 1)
        no_edges_network_after = no_edges_network - 3
        if no_edges_network > 3:
            p = (
                (move_type_probabilities[MoveType.VPLU] / move_type_probabilities[MoveType.VMIN])
                * no_edges_network
                / (no_edges_network_after**2)
            )
    if symmetries:
        # correct for number of representations, i.e., symmetries.
        result_network = result_network or apply_move(network, move)
        p *= count_automorphisms(network) / count_automorphisms(result_network)
    return p


def sample_mcmc_networks(
    number_of_leaves,
    max_retics,
    move_type_probabilities,
    correct_symmetries=True,
    burn_in=1000,
    number_of_networks=1,
):
    network = nx.DiGraph()
    network.add_nodes_from([str(k) for k in range(2 * number_of_leaves)])
    network.add_edges_from([(str(k), str(k + number_of_leaves)) for k in range(number_of_leaves)])
    current_reticulation_number = 0
    root = Root(network)
    network.add_edges_from([(0, root)])
    print(network.nodes())
    available_reticulations = ["r" + str(k + 1) for k in range(max_retics)]
    available_tree_nodes = [str(k + 2 * number_of_leaves) for k in range(max_retics)]

    for index in range(number_of_networks):
        non_moves = 0
        for j in range(burn_in):
            move = Move.random_move(
                network,
                available_tree_nodes=available_tree_nodes,
                available_reticulations=available_reticulations,
                move_type_probabilities=move_type_probabilities,
            )
            result = False
            if random.random() < AcceptanceProb(
                move,
                number_of_leaves,
                current_reticulation_number,
                max_retics,
                move_type_probabilities,
                symmetries=correct_symmetries,
                network=network,
                result_network=result,
            ):
                if not (move.move_type == MoveType.VPLU and current_reticulation_number == max_retics):
                    result = DoMove(network, move)
            if result:
                if move.move_type == [MoveType.TAIL, MoveType.HEAD]:
                    network = result
                if move.move_type == MoveType.VPLU:
                    current_reticulation_number += 1
                    available_tree_nodes.remove(move[3])
                    available_reticulations.remove(move[4])
                    network = result
                if move.move_type == MoveType.VMIN:
                    current_reticulation_number -= 1
                    available_tree_nodes += [move[1][0]]
                    available_reticulations += [move[1][1]]
                    network = result
            else:
                non_moves += 1
