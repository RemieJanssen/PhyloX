from copy import deepcopy

import numpy as np

from phylox.base import find_unused_node, suppress_node
from phylox.exceptions import InvalidMoveDefinitionException, InvalidMoveException
from phylox.rearrangement.invertsequence import from_edge
from phylox.rearrangement.movability import check_valid
from phylox.rearrangement.movetype import MoveType


def apply_move(network, move):
    """
    Apply a move to the network, not in place.
    returns True if successful, and False otherwise.

    :param network: a phylogenetic network (phylox.DiNetwork).
    :param move: a move (phylox.rearrangement.move.Move) to apply to the network.
    :return: a new phylogenetic network (phylox.DiNetwork) with the move applied.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.move import apply_move, Move
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ... )
    >>> move = Move(
    ...     move_type=MoveType.HEAD,
    ...     origin=(2,5),
    ...     moving_edge=(1,3),
    ...     target=(2,4),
    ... )
    >>> new_network = apply_move(network, move)
    >>> edges = set(new_network.edges())
    >>> edges == {(0, 1), (1, 2), (1, 3), (2, 3), (3, 4), (2, 5)}
    True
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
        if hasattr(new_network, "_reticulation_number"):
            new_network._reticulation_number += 1
        return new_network
    elif move.move_type in [MoveType.VMIN]:
        new_network.remove_edge(*move.removed_edge)
        suppress_node(new_network, move.removed_edge[0])
        suppress_node(new_network, move.removed_edge[1])
        if hasattr(new_network, "_reticulation_number"):
            new_network._reticulation_number -= 1
        return new_network
    elif move.move_type in [MoveType.NONE]:
        return network
    raise InvalidMoveException("only tail or head moves are currently valid.")


def apply_move_sequence(network, seq_moves):
    """
    Apply a sequence of moves to the network, not in place.
    returns the resulting network.

    :param network: a phylogenetic network (phylox.DiNetwork).
    :param seq_moves: a sequence of moves (list of phylox.rearrangement.move.Move) to apply to the network.
    :return: a new phylogenetic network (phylox.DiNetwork) with the moves applied.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.move import apply_move_sequence, Move
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,4),(2,5),(3,6),(3,7)],
    ... )
    >>> seq_moves = [
    ...     Move(
    ...         move_type=MoveType.TAIL,
    ...         moving_edge=(3,6),
    ...         target=(0,1),
    ...         origin=(1,7),
    ...     ),
    ...     Move(
    ...         move_type=MoveType.VPLU,
    ...         start_edge=(1,2),
    ...         end_edge=(1,7),
    ...         network=network,
    ...     ),
    ... ]
    >>> new_network = apply_move_sequence(network, seq_moves)
    >>> edges = set(new_network.edges())
    >>> edges == {(0, 3), (3, 1), (1, -1), (-1, 2), (1, -2), (-2, 7), (2, 4), (2, 5), (3, 6), (-1, -2)}
    True
    """

    for move in seq_moves:
        network = apply_move(network, move)
    return network


class Move(object):
    """
    A move is a rearrangement operation on a phylogenetic network.

    :param move_type: the type of move (phylox.rearrangement.movetype.MoveType).
    :param kwargs: the parameters of the move, depending on the move type.
    """

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
        """
        Checks if the move is of a given type.

        :param move_type: the move type to check (phylox.rearrangement.movetype.MoveType).
        :return: True if the move is of the given type, False otherwise.

        :example:
        >>> from phylox.rearrangement.move import Move
        >>> move = Move(
        ...     move_type=MoveType.TAIL,
        ...     origin=(2,5),
        ...     moving_edge=(1,3),
        ...     target=(2,4),
        ... )
        >>> move.is_type(MoveType.TAIL)
        True
        >>> move.is_type(MoveType.RSPR)
        True
        >>> move.is_type(MoveType.HEAD)
        False
        >>> move.is_type(MoveType.VERT)
        False
        """

        if move_type == MoveType.ALL:
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
        """
        Generates a random move for the given network.
        The move may not be valid.

        :param network: a phylogenetic network (phylox.DiNetwork).
        :param available_tree_nodes: a list of available tree nodes to use for the move.
        :param available_reticulations: a list of available reticulations to use for the move.
        :param move_type_probabilities: a dictionary of move type probabilities.
        :return: a random move (phylox.rearrangement.move.Move).

        :example:
        >>> from phylox import DiNetwork
        >>> from phylox.rearrangement.move import Move
        >>> network = DiNetwork(
        ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
        ... )
        >>> move = Move.random_move(
        ...     network,
        ...     move_type_probabilities={
        ...         MoveType.TAIL: 0.5,
        ...         MoveType.HEAD: 0.5,
        ...     },
        ... )
        >>> move.is_type(MoveType.RSPR)
        True
        """

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
            moving_edge_index = np.random.choice(num_edges)
            target_index = np.random.choice(num_edges - 1)
            if target_index >= moving_edge_index:
                target_index += 1
                target_index %= num_edges
            moving_edge = edges[moving_edge_index]
            target = edges[target_index]
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
