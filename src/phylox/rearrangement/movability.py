# These functions are implementations of the algorithms of Remie Janssen's PhD thesis

import networkx as nx
import random
import sys
from copy import deepcopy
from collections import deque
import re
import ast
import time

from phylox.rearrangement.move import Move
from phylox.rearrangement.movetype import MoveType


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
            raise InvalidMove(
                "origin does not match parent and child or moving_endpoint"
            )
        if network.has_edge(move.origin[0], move.origin[1]):
            raise InvalidMove("removal creates parallel edges")

        if move.is_type(MoveType.TAIL):
            if nx.has_path(network, move.moving_edge[1], move.target[0]):
                raise InvalidMove("reattachment would create a cycle")
            if move.target[1] == move.moving_edge[1]:
                raise InvalidMove("reattachment creates parallel edges")
            return
        # move.is_type(MoveType.HEAD)
        if nx.has_path(network, move.target[1], move.moving_edge[0]):
            raise InvalidMove("reattachment would create a cycle")
        if move.target[0] == move.moving_edge[0]:
            raise InvalidMove("reattachment creates parallel edges")
        return

    raise InvalidMove("Only rSPR moves are supported currently")



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
########                  Finding nodes in a network               #############
########                                                           #############
################################################################################
################################################################################
################################################################################

# Returns all nodes below a given node (including the node itself)
def AllBelow(network, node):
    """
    Finds all nodes below a given node in a network.

    :param network: a phylogenetic network.
    :param node: a node in the network.
    :return: all nodes in the network that are below the chosen node, including this node.
    """
    lengths = nx.single_source_shortest_path_length(network, node)
    return lengths.keys()


# Find the lowest nodes above a given set.
# The excluded set must include all leaves.
def LowestReticAndTreeNodeAbove(network, excludedSet, allnodes=False):
    """
    Finds a list of lowest tree nodes and a list of lowest reticulation nodes above a given set of nodes.

    :param network: a phylogenetic network.
    :param excludedSet: a set of nodes of the network.
    :param allnodes: a boolean value that determines whether we try to find all lowest nodes (True) or only one lowest node of each type (False, Default).
    :return: A list of tree nodes and a list of reticulation nodes, so that each element of these lists has all their children in the excludedSet. If not allNodes, then these lists have length at most 1.
    """
    retic = None
    tree_node = None
    lowest_retics = []
    lowest_tree_nodes = []
    for node in network.nodes():
        if node not in excludedSet:
            for c in network.successors(node):
                if c not in excludedSet:
                    break
            # else runs if the loop was not ended by a break
            # this happens exactly when all of the children are in excludedSet
            else:
                if network.out_degree(node) == 2:
                    if allnodes:
                        lowest_tree_nodes += [node]
                    elif tree_node == None:
                        # For simplicity in description, take the FIRST lowest tree node that we encounter (sim. for the reticulations)
                        tree_node = node
                elif network.in_degree(node) == 2:
                    if allnodes:
                        lowest_retics += [node]
                    elif retic == None:
                        retic = node
                if not allnodes and tree_node != None and retic != None:
                    # stop if both types of lowest nodes are found
                    break
    if allnodes:
        return lowest_tree_nodes, lowest_retics
    return tree_node, retic


# Find the highest nodes below a given set.
# The excluded set must include the root.
def HighestNodesBelow(network, excludedSet, allnodes=False):
    """
    Finds a list of highest tree nodes and a list of highest reticulation nodes below a given set of nodes.

    :param network: a phylogenetic network.
    :param excludedSet: a set of nodes of the network.
    :param allnodes: a boolean value that determines whether we try to find all highest nodes (True) or only one highest node of each type (False, Default).
    :return: A list of tree nodes and a list of reticulation nodes, so that each element of these lists has all their parents in the excludedSet. If not allNodes, then these lists have length at most 1.
    """
    retic = None
    tree_node = None
    highest_retics = []
    highest_tree_nodes = []
    for node in network.nodes():
        if node not in excludedSet:
            for c in network.predecessors(node):
                if c not in excludedSet:
                    break
            # else runs if the loop was not ended by a break
            # this happens exactly when all of the parents are in excludedSet
            else:
                if network.out_degree(node) == 2:
                    if allnodes:
                        highest_tree_nodes += [node]
                    elif tree_node == None:
                        # For simplicity in description, take the FIRST highest tree node that we encounter (sim. for the reticulations and leaves)
                        tree_node = node
                elif network.in_degree(node) == 2:
                    if allnodes:
                        highest_retics += [node]
                    elif retic == None:
                        retic = node
                if not allnodes and retic != None and tree_node != None:
                    # stop if all types of highest nodes are found
                    break
    if allnodes:
        return highest_tree_nodes, highest_retics
    return tree_node, retic


def FindTreeNode(network, excludedSet=[], randomNodes=False):
    """
    Finds a (random) tree node in a network.

    :param network: a phylogenetic network.
    :param excludedSet: a set of nodes of the network.
    :param randomNodes: a boolean value.
    :return: a tree node of the network not in the excludedSet, or None if no such node exists. If randomNodes, then a tree node is selected from all candidates uniformly at random.
    """
    all_found = []
    for node in network.nodes():
        if node not in excludedSet and network.out_degree(node) == 2 and network.in_degree(node) == 1:
            if randomNodes:
                all_found += [node]
            else:
                return node
    if all_found and randomNodes:
        return random.choice(all_found)
    return None


def FindLeaf(network, excludedSet=[], excludedParents=[], randomNodes=False):
    """
    Finds a (random) leaf in a network.

    :param network: a phylogenetic network.
    :param excludedSet: a set of nodes of the network.
    :param excludedParents: a set of nodes of the network.
    :param randomNodes: a boolean value.
    :return: a leaf of the network not in the excludedSet so that its parent is nt in excludedParents, or None if no such node exists. If randomNodes, then a leaf is selected from all candidates uniformly at random.
    """
    all_found = []
    for node in network.nodes():
        parent = Parent(network, node)
        if network.out_degree(node) == 0 and parent not in excludedParents and node not in excludedSet:
            if randomNodes:
                all_found += [node]
            else:
                return node
    if all_found and randomNodes:
        return random.choice(all_found)
    return None


def FindRetic(network, excludedSet=[], randomNodes=False):
    """
    Finds a (random) reticulation in a network.

    :param network: a phylogenetic network.
    :param excludedSet: a set of nodes of the network.
    :param randomNodes: a boolean value.
    :return: a reticulation node of the network not in the excludedSet, or None if no such node exists. If randomNodes, then a reticulation is selected from all candidates uniformly at random.
    """
    all_found = []
    for node in network.nodes():
        if node not in excludedSet and network.in_degree(node) == 2:
            if randomNodes:
                all_found += [node]
            else:
                return node
    if all_found and randomNodes:
        return random.choice(all_found)
    return None


def Parent(network, node, exclude=[], randomNodes=False):
    """
    Finds a parent of a node in a network.

    :param network: a phylogenetic network.
    :param node: a node in the network.
    :param exclude: a set of nodes of the network.
    :param randomNodes: a boolean value.
    :return: a parent of node that is not in the set of nodes exclude. If randomNodes, then this parent is selected uniformly at random from all candidates.
    """
    parent = None
    for p in network.predecessors(node):
        if p not in exclude:
            if not randomNodes:
                return p
            elif parent == None or random.getrandbits(1):
                # As there are at most two parents, we can simply replace the previous parent with probability .5 to get a random parent
                parent = p
    return parent


def Child(network, node, exclude=[], randomNodes=False):
    """
    Finds a child node of a node in a network.

    :param network: a phylogenetic network.
    :param node: a node in the network.
    :param exclude: a set of nodes of the network.
    :param randomNodes: a boolean value.
    :return: a child of node that is not in the set of nodes exclude. If randomNodes, then this child node is selected uniformly at random from all candidates.
    """
    child = None
    for c in network.successors(node):
        if c not in exclude:
            if not randomNodes:
                return c
            elif child == None or random.getrandbits(1):
                # As there are at most two children, we can simply replace the previous child with probability .5 to get a random parent
                child = c
    return child


# Returns the root of a network
def Root(network):
    """
    Finds the root of a phylogenetic network.

    :param network: a phylogenetic network.
    :return: the root node of this network.
    """
    for node in network.nodes():
        if network.in_degree(node) == 0:
            return node
    return None


# Returns a dictionary with node labels, keyed by the labels
def Labels(network):
    """
    Returns the correspondence between the leaves and the leaf-labels of a given network

    :param network: a phylogenetic network
    :return: a dictionary, where the keys are labels and the values are nodes of the network.
    """
    label_dict = dict()
    for node in network.nodes():
        node_label = network.node[node].get('label')
        if node_label:
            label_dict[node_label] = node
    return label_dict


################################################################################
################################################################################
################################################################################
########                                                           #############
########                 Sequence finding Functions                #############
########                                                           #############
################################################################################
################################################################################
################################################################################








# Find a sequence by choosing the move that most decreases the upper bound on the number of moves
# This works as long as we can always decrease the bound.
# E.g.1, this upper bound can be the length of the sequence given by Green_Line(N1,N2), the bound can always decrease after one move, if we take the move from the GL sequence (IMPLEMENTED)
# TODO E.g.2, take the upper bound given by this algorithm with bound Green_Line
def Deep_Dive_Scored(network1, network2, head_moves=True, bound_heuristic=Green_Line):
    """
    An experimental method that returns a sequence of tail/rSPR moves from network1 to network2, using the isomorphism-building heuristic for the chosen type of moves.


    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param head_moves: a boolean value that determines whether head moves are used in addition to tail moves. If True we use rSPR moves, if False we use only tail moves.
    :param bound_heuristic: a heuristic that finds a sequence between the networks quickly.
    :return: a sequence of moves from network1 to network2.
    """
    if Isomorphic(network1, network2):
        return []
    seq = []
    current_network = network1
    current_best = []
    for move in bound_heuristic(network1, network2, head_moves=head_moves):
        current_best += [(move[0], move[1], move[3])]
    if not current_best:
        return False
    done = False
    current_length = 0
    while not done:
        candidate_moves = AllValidMoves(current_network, tail_moves=True, head_moves=head_moves)
        for move in candidate_moves:
            candidate_network = DoMove(current_network, *move)
            if Isomorphic(candidate_network, network2):
                return current_best[:current_length] + [move]
            candidate_sequence = bound_heuristic(candidate_network, network2, head_moves=head_moves)
            if current_length + len(candidate_sequence) + 1 < len(current_best):
                current_best = current_best[:current_length] + [move]
                for move2 in candidate_sequence:
                    current_best += [(move2[0], move2[1], move2[3])]
        next_move = current_best[current_length]
        current_network = DoMove(current_network, *next_move)
        current_length += 1
    return True


def Depth_First(network1, network2, tail_moves=True, head_moves=True, max_time=False, show_bounds=True):
    """
    An implementation of Algorithm 1. Uses an iterated Depth First Search to simulate a Breath First Search.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param tail_moves: a boolean value determining whether tail moves are used.
    :param head_moves: a boolean value determining whether head moves are used.
    :param max_time: a float, a time limit for the function in seconds. If False, no time limit is used, and the function continues until it finds a sequence.
    :param show_bounds: a boolean parameter, if True the current lower bounds are printed to the terminal, used for debugging.
    :return: a shortest sequence of moves between the networks if it is found within the time limit, otherwise it returns an integer: a lower bound for the length of teh shortest sequence between the networks.
    """
    done = False
    lower_bound = 0
    stop_time = False
    if max_time:
        stop_time = time.time() + max_time
    while not done:
        output = Depth_First_Bounded(network1, network2, tail_moves=tail_moves, head_moves=head_moves,
                                     max_depth=lower_bound, stop_time=stop_time)
        if output == "timeout":
            return lower_bound
        elif type(output) == list:
            return output
        lower_bound += 1
        if show_bounds:
            print(lower_bound)


# Finds a shortest sequence between network1 and network2 using DFS with bounded depth
def Depth_First_Bounded(network1, network2, tail_moves=True, head_moves=True, max_depth=0, stop_time=False):
    """
    An subroutine of Algorithm 1. A depth-bounded Depth First Search used to simulate a Breath First Search.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param tail_moves: a boolean value determining whether tail moves are used.
    :param head_moves: a boolean value determining whether head moves are used.
    :param max_depth: a integer, the maximum depth for the search tree.
    :param stop_time: a float, a time limit for the function in clock time. If False, no time limit is used, and the function continues until it finds a sequence.
    :return: a shortest sequence of at most max_depth moves between the networks if it is found before the stop_time, otherwise it returns an False.
    """
    # If we cannot do any moves:
    if not tail_moves and not head_moves:
        if Isomorphic(network1, network2):
            return 0
        else:
            return False
    # Else, make a stack and search
    stack = [[]]
    while stack:
        current_moves = stack.pop()
        current_length = len(current_moves)
        current_network = network1
        for move in current_moves:
            current_network = DoMove(current_network, *move)
        if current_length == max_depth and Isomorphic(current_network, network2):
            return current_moves
        if current_length < max_depth:
            validMoves = AllValidMoves(current_network, tail_moves=tail_moves, head_moves=head_moves)
            for move in validMoves:
                stack.append(current_moves + [move])
        if stop_time and time.time() > stop_time:
            return "timeout"
    return False


# Finds a shortest sequence between network1 and network2 using BFS
def Breadth_First(network1, network2, tail_moves=True, head_moves=True, max_time=False):
    """
    A true BFS implementation to find a shortest sequence between two networks. This implementation uses too much memory, use Depth_First!

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param tail_moves: a boolean value determining whether tail moves are used.
    :param head_moves: a boolean value determining whether head moves are used.
    :param max_time: a float, a time limit for the function in seconds. If False, no time limit is used, and the function continues until it finds a sequence.
    :return: a shortest sequence of moves between the networks if it is found within the time limit, otherwise it returns an integer: a lower bound for the length of teh shortest sequence between the networks.
    """
    # If we cannot do any moves:
    if not tail_moves and not head_moves:
        if Isomorphic(network1, network2):
            return 0
        else:
            return False
    # Else, make a queue and search
    queue = deque([[]])
    start_time = time.time()
    while queue:
        current_moves = queue.popleft()
        current_network = network1
        for move in current_moves:
            current_network = DoMove(current_network, *move)
        if Isomorphic(current_network, network2):
            return current_moves
        validMoves = AllValidMoves(current_network, tail_moves=tail_moves, head_moves=head_moves)
        for move in validMoves:
            queue.append(current_moves + [move])
        if max_time and time.time() - start_time > max_time:
            return len(current_moves)
    return False


################################################################################
################################################################################
################################################################################
########                                                           #############
########                         Isomorphism                       #############
########                                                           #############
################################################################################
################################################################################
################################################################################


# Checks whether the nodes with the given attributes have the same label
def SameLabels(node1_attributes, node2_attributes):
    """
    Checks whether two nodes have the same label

    :param node1_attributes: the attributes of a node
    :param node2_attributes: the attributes of a node
    :return: True if the label attribute is the same, False otherwise.
    """
    return node1_attributes.get('label') == node2_attributes.get('label')


# Checks whether two networks are labeled isomorpgic
def Isomorphic(network1, network2):
    """
    Determines whether two networks are labeled isomorphic.

    :param network1: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute `label'.
    :param network2: a phylogenetic network, i.e., a DAG with leaf labels stored as the node attribute `label'.
    :return: True if the networks are labeled isomorphic, False otherwise.
    """
    return nx.is_isomorphic(network1, network2, node_match=SameLabels)


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
