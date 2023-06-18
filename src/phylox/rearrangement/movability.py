# These functions are implementations of the algorithms of Remie Janssen's PhD thesis

import networkx as nx
import random
import sys
from copy import deepcopy
from collections import deque
import re
import ast
import time

from phylox.rearrangement.movetype import MoveType
from phylox.exceptions import InvalidMoveException, InvalidMoveDefinitionException



def check_valid(network, move):
    """
    Checks whether a move is valid.

    :param move: a move of type Move
    :return: void
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
            raise InvalidMoveException(
                "Start node must not be in the network."
            )
        if move.end_node in network.nodes:
            raise InvalidMoveException(
                "End node must not be in the network."
            )
        if nx.has_path(network, move.end_edge[1], move.start_edge[0]) or move.start_edge == move.end_edge:
            raise InvalidMoveException("end node is reachable from start node")
    elif move.is_type(MoveType.VMIN):
        parent_0 = network.parent(move.removed_edge[0], exclude=[move.removed_edge[1]])
        child_0 = network.child(move.removed_edge[0], exclude=[move.removed_edge[1]])
        parent_1 = network.parent(move.removed_edge[1], exclude=[move.removed_edge[0]])
        child_1 = network.child(move.removed_edge[1], exclude=[move.removed_edge[0]])
        if parent_0==parent_1 and child_0==child_1:
            raise InvalidMoveException("removal creates parallel edges")
        if not (check_movable(network, move.removed_edge, move.removed_edge[0]) and check_movable(network, move.removed_edge, move.removed_edge[1])):
            raise InvalidMoveException("removal creates parallel edges")
        if network.out_degree(move.removed_edge[1]) == 0:
            raise InvalidMoveException("removes a leaf")
    else:
        raise InvalidMoveException("Only rSPR and vertical moves are supported currently")



# Returns all valid moves in the network
# List of moves in format (moving_edge,moving_endpoint,to_edge)
# change move options to move_type
def all_valid_moves(network, tail_moves=True, head_moves=True):
    """
    Finds all valid moves of a certain type in a network.

    :param network: a phylogenetic network.
    :param tail_moves: boolean value that determines whether tail moves are used (Default: True).
    :param head_moves: boolean value that determines whether head moves are used (Default: True).
    :return: A list of all valid tail/head/rSPR moves in the network.
    """
    validMoves = []
    headOrTail = []
    if tail_moves:
        headOrTail.append(0)
    if head_moves:
        headOrTail.append(1)
    for moving_edge in network.edges():
        for to_edge in network.edges():
            for end in headOrTail:
                if CheckValid(network, moving_edge, moving_edge[end], to_edge) and moving_edge[
                    end] not in to_edge:  # Last part is to prevent valid moves that result in isomorphic networks. TODO: To use this for internally labeled networks, remove this condition
                    validMoves.append((moving_edge, moving_edge[end], to_edge))
    return validMoves




# Checks whether an endpoint of an edge is movable.
def check_movable(network, moving_edge, moving_endpoint):
    """
    Checks whether an endpoint of an edge is movable in a network.

    :param network: a phylogenetic network, i.e., a DAG with labeled leaves.
    :param moving_edge: an edge in the network.
    :param moving_endpoint: a node, specifically, an endpoint of the moving_edge.
    :return: True if the endpoint of the edge is movable in the network, False otherwise.
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
    parent_of_moving_endpoint = network.parent(moving_endpoint, exclude=[moving_edge[0]])
    child_of_moving_endpoint = network.child(moving_endpoint, exclude=[moving_edge[1]])
    # if there is an edge from the parent to the child, there is a triangle
    # Otherwise, it is a movable edge
    return not network.has_edge(parent_of_moving_endpoint, child_of_moving_endpoint)
