# These functions are implementations of the algorithms of Remie Janssen's PhD thesis

import ast
import random
import re
import sys
import time
from collections import deque
from copy import deepcopy

import networkx as nx

from phylox.exceptions import InvalidMoveDefinitionException, InvalidMoveException
from phylox.rearrangement.movetype import MoveType


def check_valid(network, move):
    """
    Checks whether a move is valid.

    :param move: a rearrangement move (see phylox.rearrangement.movetype.Move)
    :return: void
    :exception: InvalidMoveException if the move is not valid

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.move import Move
    >>> from phylox.rearrangement.movability import check_valid
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ... )
    >>> move = Move(
    ...     move_type=MoveType.HEAD,
    ...     origin=(2, 5),
    ...     moving_edge=(1, 3),
    ...     target=(2, 4),
    ... )
    >>> check_valid(network, move)
    """
    if move.move_type == MoveType.NONE:
        return
    if move.is_type(MoveType.RSPR):
        if not network.has_edge(
            move.origin[0], move.moving_node
        ) or not network.has_edge(move.moving_node, move.origin[1]):
            # also catches wrong endpoint type:
            # e.g.: reticulation moving_node for tail move
            raise InvalidMoveException(
                "origin does not match parent and child or moving_endpoint"
            )
        if network.has_edge(move.origin[0], move.origin[1]):
            raise InvalidMoveException("removal creates parallel edges")

        if move.is_type(MoveType.TAIL):
            if nx.has_path(network, move.moving_edge[1], move.target[0]):
                raise InvalidMoveException("reattachment would create a cycle")
            if move.target[1] == move.moving_edge[1]:
                raise InvalidMoveException("reattachment creates parallel edges")
            return
        # move.is_type(MoveType.HEAD)
        if nx.has_path(network, move.target[1], move.moving_edge[0]):
            raise InvalidMoveException("reattachment would create a cycle")
        if move.target[0] == move.moving_edge[0]:
            raise InvalidMoveException("reattachment creates parallel edges")
    elif move.is_type(MoveType.VPLU):
        if move.start_node in network.nodes:
            raise InvalidMoveException("Start node must not be in the network.")
        if move.end_node in network.nodes:
            raise InvalidMoveException("End node must not be in the network.")
        if (
            nx.has_path(network, move.end_edge[1], move.start_edge[0])
            or move.start_edge == move.end_edge
        ):
            raise InvalidMoveException("end node is reachable from start node")
    elif move.is_type(MoveType.VMIN):
        parent_0 = network.parent(move.removed_edge[0], exclude=[move.removed_edge[1]])
        child_0 = network.child(move.removed_edge[0], exclude=[move.removed_edge[1]])
        parent_1 = network.parent(move.removed_edge[1], exclude=[move.removed_edge[0]])
        child_1 = network.child(move.removed_edge[1], exclude=[move.removed_edge[0]])
        if parent_0 == parent_1 and child_0 == child_1:
            raise InvalidMoveException("removal creates parallel edges")
        if not (
            check_movable(network, move.removed_edge, move.removed_edge[0])
            and check_movable(network, move.removed_edge, move.removed_edge[1])
        ):
            raise InvalidMoveException("removal creates parallel edges")
        if network.out_degree(move.removed_edge[1]) == 0:
            raise InvalidMoveException("removes a leaf")
    else:
        raise InvalidMoveException(
            "Only rSPR and vertical moves are supported currently"
        )


# Checks whether an endpoint of an edge is movable.
def check_movable(network, moving_edge, moving_endpoint):
    """
    Checks whether an endpoint of an edge is movable in a network.

    :param network: a phylogenetic network, i.e., a DAG with labeled leaves.
    :param moving_edge: an edge in the network.
    :param moving_endpoint: a node, specifically, an endpoint of the moving_edge.
    :return: True if the endpoint of the edge is movable in the network, False otherwise.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.movability import check_movable
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ... )
    >>> check_movable(network, (1, 3), 3)
    True
    >>> check_movable(network, (1, 3), 1)
    True
    >>> check_movable(network, (3, 5), 3)
    False
    """
    if moving_endpoint == moving_edge[0]:
        # Tail move
        if network.in_degree(moving_endpoint) in (0, 2):
            # cannot move the tail if it is a reticulation or root
            return False
    elif moving_endpoint == moving_edge[1]:
        # Head move
        if network.out_degree(moving_endpoint) in (0, 2):
            # cannot move the head if it is a tree node or leaf
            return False
    else:
        # Moving endpoint is not part of the moving edge
        return False
    # Now check for triangles, by finding the other parent and child of the moving endpoint
    parent_of_moving_endpoint = network.parent(
        moving_endpoint, exclude=[moving_edge[0]]
    )
    child_of_moving_endpoint = network.child(moving_endpoint, exclude=[moving_edge[1]])
    # if there is an edge from the parent to the child, there is a triangle
    # Otherwise, it is a movable edge
    return not network.has_edge(parent_of_moving_endpoint, child_of_moving_endpoint)
