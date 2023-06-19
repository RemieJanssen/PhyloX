from copy import deepcopy

import numpy as np

from phylox.base import find_unused_node
from phylox.exceptions import InvalidMoveDefinitionException, InvalidMoveException
from phylox.rearrangement.invertsequence import from_edge
from phylox.rearrangement.movability import check_valid
from phylox.rearrangement.movetype import MoveType


def apply_move(network, move):
    """
    Apply a move to the network, not in place.
    returns True if successful, and False otherwise.
    """
    check_valid(network, move)
    new_network = deepcopy(network)

    if move.move_type in [MoveType.TAIL, MoveType.HEAD]:
        if move.moving_node in move.target:
            # move does not impact the network
            return network
        new_network.remove_edges_from(
            [
                (move.origin[0], move.moving_node),
                (move.moving_node, move.origin[1]),
                move.target,
            ]
        )
        new_network.add_edges_from(
            [
                (move.target[0], move.moving_node),
                (move.moving_node, move.target[1]),
                move.origin,
            ]
        )
        return new_network
    elif move.move_type in [MoveType.VPLU]:
        new_network.remove_edges_from([move.start_edge, move.end_edge])
        new_network.add_edges_from(
            [
                (move.start_edge[0], move.start_node),
                (move.start_node, move.start_edge[1]),
                (move.end_edge[0], move.end_node),
                (move.end_node, move.end_edge[1]),
                (move.start_node, move.end_node),
            ]
        )
        return new_network
    elif move.move_type in [MoveType.VMIN]:
        parent_0 = network.parent(move.removed_edge[0], exclude=[move.removed_edge[1]])
        child_0 = network.child(move.removed_edge[0], exclude=[move.removed_edge[1]])
        parent_1 = network.parent(move.removed_edge[1], exclude=[move.removed_edge[0]])
        child_1 = network.child(move.removed_edge[1], exclude=[move.removed_edge[0]])
        new_network.remove_edges_from(
            [
                (parent_0, move.removed_edge[0]),
                (move.removed_edge[0], child_0),
                (parent_1, move.removed_edge[1]),
                (move.removed_edge[1], child_1),
                move.removed_edge,
            ]
        )
        new_network.add_edges_from([(parent_0, child_0), (parent_1, child_1)])
        return new_network
    elif move.move_type in [MoveType.NONE]:
        return network
    raise InvalidMoveException("only tail or head moves are currently valid.")


def apply_move_sequence(network, seq_moves):
    for move in seq_moves:
        network = apply_move(network, move)
    return network


class Move(object):
    def __init__(self, *args, **kwargs):
        try:
            self.move_type = kwargs["move_type"]
        except KeyError:
            raise InvalidMoveDefinitionException("Missing move_type.")

        # None type move
        if self.move_type == MoveType.NONE:
            return

        # TAIL/HEAD move (i.e. RSPR/horizontal)
        elif self.move_type == MoveType.RSPR:
            raise InvalidMoveDefinitionException(
                "rSPR moves must be defined as moves of type tail or head."
            )
        elif self.move_type in [MoveType.TAIL, MoveType.HEAD]:
            try:
                self.origin = kwargs["origin"]
                self.moving_edge = kwargs["moving_edge"]
                self.target = kwargs["target"]
            except KeyError:
                raise InvalidMoveDefinitionException(
                    "Missing one of origin, moving_edge, or target."
                )

            if self.move_type == MoveType.TAIL:
                self.moving_node = self.moving_edge[0]
            else:
                self.moving_node = self.moving_edge[1]

            if self.moving_edge == self.target:
                raise InvalidMoveDefinitionException(
                    "Moving edge must not be the target edge."
                )

            return
        # VERT move (i.e. SPR/vertical)
        elif self.move_type == MoveType.VERT:
            raise InvalidMoveDefinitionException(
                "vertical moves must be defined as moves of type VPLU or VMIN."
            )
        # VPLU/VMIN move
        elif self.move_type == MoveType.VPLU:
            try:
                self.start_edge = kwargs["start_edge"]
                self.end_edge = kwargs["end_edge"]
                self.start_node = kwargs.get("start_node", None)
                self.end_node = kwargs.get("end_node", None)
                network = kwargs.get("network", None)
                if (
                    self.start_node is None or self.end_node is None
                ) and network is None:
                    raise InvalidMoveDefinitionException(
                        "Either a start_node and end_node, or a network must be given."
                    )
                if self.start_node is None:
                    self.start_node = find_unused_node(network)
                if self.end_node is None:
                    self.end_node = find_unused_node(network, exclude=[self.start_node])

                if self.start_edge == self.end_edge:
                    raise InvalidMoveDefinitionException(
                        "Start edge must not be the end edge."
                    )
            except KeyError:
                raise InvalidMoveDefinitionException(
                    "Missing one of start_edge, end_edge."
                )
        elif self.move_type == MoveType.VMIN:
            try:
                self.removed_edge = kwargs["removed_edge"]
            except KeyError:
                raise InvalidMoveDefinitionException(
                    "Missing removed_edge in definition."
                )
        else:
            raise InvalidMoveDefinitionException("Invalid move type.")

    def is_type(self, move_type):
        if self.move_type == MoveType.ALL:
            return True
        if (
            self.move_type == MoveType.NONE
            or (
                move_type == MoveType.RSPR
                and self.move_type in [MoveType.TAIL, MoveType.HEAD]
            )
            or (
                move_type == MoveType.VERT
                and self.move_type in [MoveType.VPLU, MoveType.VMIN]
            )
        ):
            return True
        return move_type == self.move_type

    @classmethod
    def random_move(
        cls,
        network,
        available_tree_nodes=None,
        available_reticulations=None,
        move_type_probabilities={
            MoveType.TAIL: 0.4,
            MoveType.HEAD: 0.4,
            MoveType.VPLU: 0.1,
            MoveType.VMIN: 0.1,
        },
    ):
        available_tree_nodes = available_tree_nodes or []
        available_reticulations = available_reticulations or []
        edges = list(network.edges())
        num_edges = len(edges)
        move_type_probabilities_keys = list(move_type_probabilities.keys())
        movetype_index = np.random.choice(
            len(move_type_probabilities_keys), p=tuple(move_type_probabilities.values())
        )
        movetype = move_type_probabilities_keys[movetype_index]
        if movetype in [MoveType.TAIL, MoveType.HEAD]:
            moving_edge = edges[np.random.choice(num_edges)]
            target = edges[np.random.choice(num_edges)]
            moving_endpoint_index = 0 if movetype == MoveType.TAIL else 1
            moving_endpoint = moving_edge[moving_endpoint_index]
            origin = from_edge(network, moving_edge, moving_endpoint=moving_endpoint)
            return Move(
                move_type=movetype,
                origin=origin,
                moving_edge=moving_edge,
                target=target,
            )
        elif movetype == MoveType.VPLU:
            available_reticulations = list(available_reticulations) or [
                find_unused_node(network)
            ]
            available_tree_nodes = list(available_tree_nodes) or [
                find_unused_node(network, exclude=available_reticulations)
            ]
            start_edge = edges[np.random.choice(num_edges)]
            end_edge = edges[np.random.choice(num_edges)]
            start_node = np.random.choice(available_tree_nodes)
            end_node = np.random.choice(available_reticulations)
            return Move(
                move_type=movetype,
                start_edge=start_edge,
                end_edge=end_edge,
                start_node=start_node,
                end_node=end_node,
            )
        elif movetype == MoveType.VMIN:
            removed_edge = edges[np.random.choice(num_edges)]
            return Move(move_type=movetype, removed_edge=removed_edge)
