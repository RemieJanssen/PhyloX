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
        if nx.has_path(network, move.end_edge[1], move.start_edge[0]) or move.start_edge == move.end_edge:
            raise InvalidMoveException("end node is reachable from start node")
    elif move.is_type(MoveType.VMIN):
        parent_0 = network.parent(removed_edge[0], exclude=[removed_edge[1]])
        child_0 = network.child(removed_edge[0], exclude=[removed_edge[1]])
        parent_1 = network.parent(removed_edge[1], exclude=[removed_edge[0]])
        child_1 = network.child(removed_edge[1], exclude=[removed_edge[0]])
        if parent_0==parent_1 and child_0==child_1:
            raise InvalidMoveException("removal creates parallel edges")
        if not (CheckMovable(network, move.removed_edge, move.removed_edge[0]) and CheckMovable(network, move.removed_edge, move.removed_edge[1])):
            raise InvalidMoveException("removal creates parallel edges")
    else:
        raise InvalidMoveException("Only rSPR and vertical moves are supported currently")



# Returns all valid moves in the network
# List of moves in format (moving_edge,moving_endpoint,to_edge)
def AllValidMoves(network, tail_moves=True, head_moves=True):
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






################################################################################
################################################################################
################################################################################
########                                                           #############
########                 Sequence finding Functions                #############
########                                                           #############
################################################################################
################################################################################
################################################################################














################################################################################
################################################################################
################################################################################
########                                                           #############
########                         Isomorphism                       #############
########                                                           #############
################################################################################
################################################################################
################################################################################




################################################################################
################################################################################
################################################################################
########                                                           #############
########                         I/O Functions                     #############
########                                                           #############
################################################################################
################################################################################
################################################################################


########
######## Convert Newick to a networkx Digraph with labels (and branch lengths)
########
# Write length newick: convert ":" to "," and then evaluate as list of lists using ast.literal_eval
# Then, in each list, the node is followed by the length of the incoming arc.
# This only works as long as each branch has length and all internal nodes are labeled.
def Newick_To_Network(newick):
    """
    Converts a Newick string to a networkx DAG with leaf labels.

    :param newick: a string in extended Newick format for phylogenetic networks.
    :return: a phylogenetic network, i.e., a networkx digraph with leaf labels represented by the `label' node attribute.
    """
    # Ignore the ';'
    newick = newick[:-1]
    # If names are not in string format between ', add these.
    if not "'" in newick and not '"' in newick:
        newick = re.sub(r"\)#H([\d]+)", r",#R\1)", newick)
        newick = re.sub(r"([,\(])([#a-zA-Z\d]+)", r"\1'\2", newick)
        newick = re.sub(r"([#a-zA-Z\d])([,\(\)])", r"\1'\2", newick)
    else:
        newick = re.sub(r"\)#H([d]+)", r"'#R\1'\)", newick)
    newick = newick.replace("(", "[")
    newick = newick.replace(")", "]")
    nestedtree = ast.literal_eval(newick)
    edges, current_node = NestedList_To_Network(nestedtree, 0, 1)
    network = nx.DiGraph()
    network.add_edges_from(edges)
    network = NetworkLeafToLabel(network)
    return network


# Returns a network in which the leaves have the original name as label, and all nodes have integer names.
def NetworkLeafToLabel(network):
    """
    Renames the network nodes to integers, while storing the original node names in the `original' node attribute.

    :param network: a phylogenetic network
    :return: a phylogenetic network with original node names in the `original' node attribute.
    """
    for node in network.nodes():
        if network.out_degree(node) == 0:
            network.node[node]['label'] = node
    return nx.convert_node_labels_to_integers(network, first_label=0, label_attribute='original')


# Auxiliary function to convert list of lists to graph
def NestedList_To_Network(nestedList, top_node, next_node):
    """
    Auxiliary function used to convert list of lists to graph.

    :param nestedList: a nested list.
    :param top_node: an integer, the node name of the top node of the network represented by the list.
    :param next_node: an integer, the lowest integer not yet used as node name in the network.
    :return: a set of edges of the network represented by the nested list, and an updated next_node.
    """
    edges = []
    if type(nestedList) == list:
        if type(nestedList[-1]) == str and len(nestedList[-1]) > 2 and nestedList[-1][:2] == '#R':
            retic_node = '#H' + nestedList[-1][2:]
            bottom_node = retic_node
        else:
            bottom_node = next_node
            next_node += 1
        edges.append((top_node, bottom_node))
        for t in nestedList:
            extra_edges, next_node = NestedList_To_Network(t, bottom_node, next_node)
            edges = edges + extra_edges
    else:
        if not (len(nestedList) > 2 and nestedList[:2] == '#R'):
            edges = [(top_node, nestedList)]
    return edges, next_node


# Sets the labels of the nodes of a network with a given label dictionary
def Set_Labels(network, label_dict):
    """
    Sets the labels of the leaves of a network using a dictionary of labels.

    :param network: a networkx digraph, a DAG.
    :param label_dict: a dictionary, containing the labels (values) of the nodes of the network (keys).
    :return: a phylogenetic network, obtained by labeling network with the labels.
    """
    for node, value in label_dict.items():
        network.node[node]['label'] = value
