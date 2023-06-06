# These functions are implementations of the algorithms of Remie Janssen's PhD thesis

import networkx as nx
import random
import sys
from copy import deepcopy
from collections import deque
import re
import ast
import time


# Contents
# - Move functions
# - Finding nodes in a network
# - Sequence finding Functions
# - Isomorphism
# - I/O Functions


################################################################################
################################################################################
################################################################################
########                                                           #############
########                       Move functions                      #############
########                                                           #############
################################################################################
################################################################################
################################################################################


# Checks whether an endpoint of an edge is movable.
def CheckMovable(network, moving_edge, moving_endpoint):
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
    parent_of_moving_endpoint = Parent(network, moving_endpoint, exclude=[moving_edge[0]])
    child_of_moving_endpoint = Child(network, moving_endpoint, exclude=[moving_edge[1]])
    # if there is an edge from the parent to the child, there is a triangle
    # Otherwise, it is a movable edge
    return not network.has_edge(parent_of_moving_endpoint, child_of_moving_endpoint)


# Checks whether the given move is valid.
def CheckValid(network, moving_edge, moving_endpoint, to_edge):
    """
    Checks whether a move is valid.

    :param network: a phylogenetic network, i.e., a DAG with labeled leaves.
    :param moving_edge: an edge in the network.
    :param moving_endpoint: a node, specifically, an endpoint of the moving_edge.
    :param to_edge: another edge of the network.
    :return: True if moving moving_endpoint of moving_edge to to_edge is allowed, False otherwise.
    """
    if moving_edge == to_edge:
        return False
    if not CheckMovable(network, moving_edge, moving_endpoint):
        return False
    if moving_endpoint == moving_edge[0]:
        # tail move, check whether the to_edge is below the head of the moving edge
        if nx.has_path(network, moving_edge[1], to_edge[0]):
            # the move would create a cycle
            return False
        if to_edge[1] == moving_edge[1]:
            # the move would create parallel edges
            return False
    elif moving_endpoint == moving_edge[1]:
        # head move, check whether the to_edge is above the tail of the moving edge
        if nx.has_path(network, to_edge[1], moving_edge[0]):
            # the move would create a cycle
            return False
        if to_edge[0] == moving_edge[0]:
            # the move would create parallel edges
            return False
    else:
        # The moving endpoint is not part of the moving edge
        # Checked in CheckMovable as well, redundant?!
        return False
    # No parallel edges at start location
    # No cycle
    # No parallel edges at end location
    # So the move is valid
    return True


# Attempts to do the move of given moving_edge (with a given endpoint, the head or the tail), to the edge to_edge.
# if check_validity==False, the validity of the move will not be checked, and the move will be performed whether it is possible or not.
# To use this setting, make sure that the move is valid yourself!
# Returns False if the move is invalid
# Returns the new network if the move is valid.
def DoMove(network, moving_edge, moving_endpoint, to_edge, check_validity=True):
    """
    Performs a move on a network and returns the result (returns False if validity is checked to be False).

    :param network: a phylogenetic network, i.e., a DAG with labeled leaves.
    :param moving_edge: an edge in the network.
    :param moving_endpoint: a node, specifically, an endpoint of the moving_edge.
    :param to_edge: another edge of the network.
    :param check_validity: If True, performs the move only after checking the validity of the move. If False, the move is performed even if it is invalid, which may result in a graph that is not a phylogenetic network.
    :return: False if validity is checked and the move is not valid. Otherwise, returns the resulting network of the move.
    """
    # Check wether the move is valid if the user wants to.
    if check_validity:
        if not CheckValid(network, moving_edge, moving_endpoint, to_edge):
            return False
    from_edge = From_Edge(network, moving_edge, moving_endpoint)
    # Perform the move and return the result
    # If the move is to an edge adjacent to the moving endpoint, do nothing
    if moving_endpoint in to_edge:
        return network.copy()

    # else, do the actual moves.
    new_network = network.copy()
    new_network.remove_edges_from([(from_edge[0], moving_endpoint), (moving_endpoint, from_edge[1]), to_edge])
    new_network.add_edges_from([(to_edge[0], moving_endpoint), (moving_endpoint, to_edge[1]), from_edge])
    return new_network


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


def GL_Case1_rSPR(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False):
    """
    An implementation of Algorithm 2. Finds a sequence of rSPR moves that makes it possible to add the lowest reticulation up to the down-closed isomrophism.

    :param N: a phylogenetic network.
    :param Np: a phylogenetic network.
    :param up: a lowest reticulation node of Np above the isomorphism.
    :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
    :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
    :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
    :return: a list of rSPR moves in N, a list of rSPR moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism.
    """
    # use notation as in the paper
    # ' is denoted p
    xp = Child(Np, up)
    x = isom_Np_N[xp]
    z = Parent(N, x, exclude=isom_N_Np.keys(), randomNodes=randomNodes)

    # Case1a: z is a reticulation
    if N.in_degree(z) == 2:
        return [], [], z, up
    # Case1b: z is not a reticulation
    # Find a retic v in N not in the isom yet
    v = FindRetic(N, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes)
    u = None
    for parent in N.predecessors(v):
        if CheckValid(N, (parent, v), v, (z, x)):
            if not randomNodes:
                u = parent
                break
            elif u == None or random.getrandbits(1):
                u = parent
    # If v has an incoming arc (u,v) movable to (z,x), perform this move
    if u != None:
        return [((u, v), v, (z, x))], [], v, up
    # if none of the moves is valid, this must be because
    # v is already a reticulation above x with its movable incoming arc (u,v) with u=z.
    return [], [], v, up


def GL_Case1_Tail(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False):
    """
    An implementation of Algorithm 6. Finds a sequence of tail moves that makes it possible to add the lowest reticulation up to the down-closed isomrophism.

    :param N: a phylogenetic network.
    :param Np: a phylogenetic network.
    :param up: a lowest reticulation node of Np above the isomorphism.
    :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
    :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
    :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
    :return: a list of tail moves in N, a list of tail moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism. Returns false if the networks are not isomorphic with 2 leaves and 1 reticulation.
    """
    # use notation as in the paper
    # ' is denoted p
    xp = Child(Np, up)
    x = isom_Np_N[xp]
    z = Parent(N, x, exclude=isom_N_Np.keys(), randomNodes=randomNodes)
    # Case1a: z is a reticulation
    if N.in_degree(z) == 2:
        return [], [], z, up
    # Case1b: z is not a reticulation
    # z is a tree node
    if CheckMovable(N, (z, x), z):
        # Case1bi: (z,x) is movable
        # Find a reticulation u in N not in the isomorphism yet
        # TODO: Can first check if the other parent of x suffices here, should heuristcally be better
        u = FindRetic(N, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes)
        v = Child(N, u)
        if v == x:
            return [], [], u, up
        # we may now assume v!=x
        if z == v:
            v = Child(N, z, exclude=[x])
            w = Parent(N, u, randomNodes=randomNodes)
            return [((z, v), z, (w, u))], [], u, up
        w = Parent(N, u, exclude=[z], randomNodes=randomNodes)
        return [((z, x), z, (u, v)), ((z, v), z, (w, u))], [], u, up
    # Case1bii: (z,x) is not movable
    c = Parent(N, z)
    d = Child(N, z, exclude=[x])
    # TODO: b does not have to exist if we have an outdeg-2 root, this could be c!
    b = Parent(N, c)
    if N.in_degree(b) != 0:
        # Case1biiA: b is not the root of N
        a = Parent(N, b, randomNodes=randomNodes)
        # First do the move ((c,d),c,(a,b)), then Case1bi applies as (z,x) is now movable
        newN = DoMove(N, (c, d), c, (a, b))
        u = FindRetic(newN, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes)
        v = Child(newN, u)
        if v == x:
            # In this case, u is a reticulation parent of x and u is not in the isom. Hence, we can simply add it to the isom.
            # Note: The tail move we did is not necessary!
            # TODO: First check both parents of x for being a reticulation not in the isomorphism yet
            return [], [], u, up
        # we may now assume v!=x
        if z == v:
            # This only happens if z==v and u==b
            # we move z up above the retic b as well, too
            w = Parent(newN, b, randomNodes=randomNodes)
            return [((c, d), c, (a, b)), ((z, d), z, (w, b))], [], u, up
        w = Parent(newN, u, exclude=[z], randomNodes=randomNodes)
        return [((c, d), c, (a, b)), ((z, x), z, (u, v)), ((z, v), z, (w, u))], [], u, up
    # Case1biiB: b is the root of N
    # Note: d is not in the isomorphism
    e = Child(N, d)
    if e == x:
        return [], [], d, up
    if N.out_degree(x) == 2:
        s = Child(N, x)
        t = Child(N, x, exclude=[s])
        if s == e:
            return [((x, t), x, (d, e))], [], d, up
        if t == e:
            return [((x, s), x, (d, e))], [], d, up
        return [((x, s), x, (d, e)), ((x, e), x, (z, t)), ((x, t), x, (d, s))], [], d, up
    if N.out_degree(e) == 2:
        s = Child(N, e)
        t = Child(N, e, exclude=[s])
        if s == x:
            return [((e, t), e, (z, x))], [], d, up
        if t == x:
            return [((e, s), e, (z, x))], [], d, up
        return [((e, s), e, (z, x)), ((e, x), e, (d, t)), ((e, t), e, (z, s))], [], d, up
    # neither are tree nodes, so both must be leaves
    # In that case, there is no sequence between the two networks.
    return [], [], None, None


def GL_Case3(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False):
    """
    An implementation of Algorithm 3. Finds a sequence of tail moves that makes it possible to add the lowest tree node up to the down-closed isomrophism.

    :param N: a phylogenetic network.
    :param Np: a phylogenetic network.
    :param up: a lowest tree node of Np above the isomorphism.
    :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
    :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
    :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
    :return: a list of tail moves in N, a list of tail moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism.
    """
    # Find the children x' and y' of u'
    xp, yp = list(Np.successors(up))
    # Find the corresponding nodes x and y in N
    x = isom_Np_N[xp]
    y = isom_Np_N[yp]
    # Find the set of common parents of x and y
    parents_x = set(N.predecessors(x))
    parents_y = set(N.predecessors(y))
    common_parents = parents_x & parents_y
    # Case3a: x and y have a common parent not in the isom
    common_parents_not_Y = []
    for parent in common_parents:
        if parent not in isom_N_Np.keys():
            # then parent can be mapped to up
            common_parents_not_Y += [parent]
            if not randomNodes:
                return [], [], parent, up
    if common_parents_not_Y:
        return [], [], random.choice(common_parents_not_Y), up

    # Case3b: x and y do not have a common parent in the isomorphism
    # For both, find a parent not in the isomorphism yet
    # TODO: preferably make them tree nodes
    z_x = Parent(N, x, exclude=isom_N_Np.keys(), randomNodes=randomNodes)
    z_y = Parent(N, y, exclude=isom_N_Np.keys(), randomNodes=randomNodes)

    # Case3bi: (z_x,x) is movable
    if CheckValid(N, (z_x, x), z_x, (z_y, y)):
        return [((z_x, x), z_x, (z_y, y))], [], z_x, up
    # Case3bii: (z_y,y) is movable
    if CheckValid(N, (z_y, y), z_y, (z_x, x)):
        return [((z_y, y), z_y, (z_x, x))], [], z_y, up
    # Case3biii: Neither (z_x,x) nor (z_y,y) is movable

    if N.in_degree(z_x) == 2 or N.in_degree(z_y) == 2:
        return [], [], None, None
    # Both these parents are tree nodes.
    # This happens always in the non-random version, as otherwise there is a lowest reticulation node and we would be in Case 1 or 2.

    # As both nodes are tree nodes and the arcs immovable, both arcs hang of the side of a triangle.
    # Find the top node of the triangle for z_x
    c_x = Parent(N, z_x)
    b_x = Parent(N, c_x)

    # Find the top node of the triangle for z_y
    c_y = Parent(N, z_y)
    #    print()
    #    print(z_y,N.edges())

    b_y = Parent(N, c_y)

    if N.in_degree(b_x) == 0:
        # c_x is the child of the root
        # c_x!=c_y, so c_y is not the child of the root
        # swap the roles of x and y
        x, y = y, x
        z_x, z_y = z_y, z_x
        b_x, b_y = b_y, b_x
        c_x, c_y = c_y, c_x
    # c_x is not the child of the root
    # find a parent of b_x, and the bottom node of the triangle d_x
    a_x = Parent(N, b_x, randomNodes=randomNodes)
    d_x = Child(N, c_x, exclude=[z_x])
    return [((c_x, d_x), c_x, (a_x, b_x)), ((z_x, x), z_x, (z_y, y))], [], z_x, up


def Green_Line(network1, network2, head_moves=True):
    """
    An implementation of Algorithm 4 and its tail move counterpart. Finds a sequence of tail/rSPR moves from network1 to network2 by building a down-closed isomorphism.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param head_moves: a boolean value determining whether head moves are allowed (if True we use rSPR moves, if False we only use tail moves).
    :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
    """
    # Find the root and labels of the networks
    root1 = Root(network1)
    root2 = Root(network2)
    label_dict_1 = Labels(network1)
    label_dict_2 = Labels(network2)

    # initialize isomorphism
    isom_1_2 = dict()
    isom_2_1 = dict()
    isom_size = 0
    for label, node1 in label_dict_1.items():
        node2 = label_dict_2[label]
        isom_1_2[node1] = node2
        isom_2_1[node2] = node1
        isom_size += 1

    # Keep track of the size of the isomorphism and the size it is at the end of the green line algorithm
    goal_size = len(network1) - 1

    # init lists of sequence of moves
    # list of (moving_edge,moving_endpoint,from_edge,to_edge)
    seq_from_1 = []
    seq_from_2 = []
    # TODO keep track of lowest nodes? (Currently, the code does not do this, but it could speed up the code)

    # Do the green line algorithm
    while (isom_size < goal_size):
        # Find lowest nodes above the isom in the networks:
        lowest_tree_node_network1, lowest_retic_network1 = LowestReticAndTreeNodeAbove(network1, isom_1_2.keys())
        lowest_tree_node_network2, lowest_retic_network2 = LowestReticAndTreeNodeAbove(network2, isom_2_1.keys())

        ######################################
        # Case1: a lowest retic in network1
        if lowest_retic_network1 != None:
            # use notation as in the paper network1 = N', network2 = N, where ' is denoted p
            up = lowest_retic_network1
            if head_moves:
                moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_rSPR(network2,
                                                                                                             network1,
                                                                                                             up,
                                                                                                             isom_2_1,
                                                                                                             isom_1_2)
            else:
                moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_Tail(network2,
                                                                                                             network1,
                                                                                                             up,
                                                                                                             isom_2_1,
                                                                                                             isom_1_2)
                if added_node_network_1 == None:
                    return False
        ######################################
        # Case2: a lowest retic in network2
        elif lowest_retic_network2 != None:
            # use notation as in the paper network2 = N', network1 = N, where ' is denoted p
            up = lowest_retic_network2
            if head_moves:
                moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_rSPR(network1,
                                                                                                             network2,
                                                                                                             up,
                                                                                                             isom_1_2,
                                                                                                             isom_2_1)
            else:
                moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_Tail(network1,
                                                                                                             network2,
                                                                                                             up,
                                                                                                             isom_1_2,
                                                                                                             isom_2_1)
                if added_node_network_1 == None:
                    return False

                    ######################################
        # Case3: a lowest tree node in network1
        else:
            # use notation as in the paper network1 = N, network2 = N'
            up = lowest_tree_node_network2
            moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case3(network1, network2,
                                                                                                    up, isom_1_2,
                                                                                                    isom_2_1)
        # Now perform the moves and update the isomorphism
        isom_1_2[added_node_network_1] = added_node_network_2
        isom_2_1[added_node_network_2] = added_node_network_1
        for move in moves_network_1:
            seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
            network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
        for move in moves_network_2:
            seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
            network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
        isom_size += 1
    # TESTING FOR CORRECTNESS WHILE RUNNING
    #        if not Isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
    #            print("not unlabeled isom")
    #            print(seq_from_1)
    #            print(seq_from_2)
    #            print(network1.subgraph(isom_1_2.keys()).edges())
    #            print(network2.subgraph(isom_2_1.keys()).edges())

    # Add the root to the isomorphism, if it was there
    isom_1_2[root1] = root2
    isom_2_1[root2] = root1
    # invert seq_from_2, rename to node names of network1, and append to seq_from_1
    return seq_from_1 + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)


def Green_Line_Random(network1, network2, head_moves=True, repeats=1):
    """
    Finds a sequence of tail/rSPR moves from network1 to network2 by randomly building a down-closed isomorphism a number of times, and only keeping the shortest sequence.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param head_moves: a boolean value determining whether head moves are allowed (if True we use rSPR moves, if False we only use tail moves).
    :param repeats: an integer, determining how many repeats of Green_Line_Random_Single are performed.
    :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
    """
    best_seq = None
    for i in range(repeats):
        candidate_seq = Green_Line_Random_Single(network1, network2, head_moves=head_moves)
        if candidate_seq == False:
            return False
        if best_seq == None or len(best_seq) > len(candidate_seq):
            best_seq = candidate_seq
    return best_seq


def Green_Line_Random_Single(network1, network2, head_moves=True):
    """
    An implementation of Algorithm 5 and its tail move counterpart. Finds a sequence of tail/rSPR moves from network1 to network2 by randomly building a down-closed isomorphism.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param head_moves: a boolean value determining whether head moves are allowed (if True we use rSPR moves, if False we only use tail moves).
    :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
    """
    # Find the root and labels of the networks
    root1 = Root(network1)
    root2 = Root(network2)
    label_dict_1 = Labels(network1)
    label_dict_2 = Labels(network2)

    # initialize isomorphism
    isom_1_2 = dict()
    isom_2_1 = dict()
    isom_size = 0
    for label, node1 in label_dict_1.items():
        node2 = label_dict_2[label]
        isom_1_2[node1] = node2
        isom_2_1[node2] = node1
        isom_size += 1

    # Keep track of the size of the isomorphism and the size it is at the end of the green line algorithm
    goal_size = len(network1) - 1

    # init lists of sequence of moves
    # list of (moving_edge,moving_endpoint,from_edge,to_edge)
    seq_from_1 = []
    seq_from_2 = []
    # TODO keep track of lowest nodes?

    # Do the green line algorithm
    while (isom_size < goal_size):
        # Find all lowest nodes above the isom in the networks:
        lowest_tree_node_network1, lowest_retic_network1 = LowestReticAndTreeNodeAbove(network1, isom_1_2.keys(),
                                                                                       allnodes=True)
        lowest_tree_node_network2, lowest_retic_network2 = LowestReticAndTreeNodeAbove(network2, isom_2_1.keys(),
                                                                                       allnodes=True)

        # Construct a list of all lowest nodes in a tuple with the corresponding network (in random order)
        # I.e. If u is a lowest node of network one, it will appear in the list as (u,1)
        lowest_nodes_network1 = [(u, 1) for u in lowest_tree_node_network1 + lowest_retic_network1]
        lowest_nodes_network2 = [(u, 2) for u in lowest_tree_node_network2 + lowest_retic_network2]
        candidate_lowest_nodes = lowest_nodes_network1 + lowest_nodes_network2
        random.shuffle(candidate_lowest_nodes)

        # As some cases do not give an addition to the isom, we continue trying lowest nodes until we find one that does.
        for lowest_node, network_number in candidate_lowest_nodes:
            ######################################
            # Case1: a lowest retic in network1
            if network_number == 1 and network1.in_degree(lowest_node) == 2:
                # use notation as in the paper network1 = N', network2 = N, where ' is denoted p
                up = lowest_node
                if head_moves:
                    moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_rSPR(
                        network2, network1, up, isom_2_1, isom_1_2, randomNodes=True)
                else:
                    moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_Tail(
                        network2, network1, up, isom_2_1, isom_1_2, randomNodes=True)
                    if added_node_network_1 == None:
                        # The networks are non-isom networks with 2 leaves and 1 reticulation
                        return False
                # This case always gives a node to add to the isom
                break

            ######################################
            # Case2: a lowest retic in network2
            elif network_number == 2 and network2.in_degree(lowest_node) == 2:
                # use notation as in the paper network2 = N', network1 = N, where ' is denoted p
                up = lowest_node
                if head_moves:
                    moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_rSPR(
                        network1, network2, up, isom_1_2, isom_2_1, randomNodes=True)
                else:
                    moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_Tail(
                        network1, network2, up, isom_1_2, isom_2_1, randomNodes=True)
                    if added_node_network_1 == None:
                        # The networks are non-isom networks with 2 leaves and 1 reticulation
                        return False
                        # This case always gives a node to add to the isom
                break

            ######################################
            # Case3: a lowest tree node in network1
            elif network_number == 2 and network2.out_degree(lowest_node) == 2:
                # use notation as in the paper network1 = N, network2 = N'
                up = lowest_node
                moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case3(network1,
                                                                                                        network2, up,
                                                                                                        isom_1_2,
                                                                                                        isom_2_1,
                                                                                                        randomNodes=True)
                # If we can add a node to the isom, added_node_network_2 has a value
                if added_node_network_2:
                    break

            ######################################
            # Case3': a lowest tree node in network2
            else:
                # use notation as in the paper network1 = N, network2 = N'
                up = lowest_node
                moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case3(network2,
                                                                                                        network1, up,
                                                                                                        isom_2_1,
                                                                                                        isom_1_2,
                                                                                                        randomNodes=True)
                # If we can add a node to the isom, added_node_network_2 has a value
                if added_node_network_2:
                    break

        # Now perform the moves and update the isomorphism
        isom_1_2[added_node_network_1] = added_node_network_2
        isom_2_1[added_node_network_2] = added_node_network_1
        for move in moves_network_1:
            seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
            network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
        for move in moves_network_2:
            seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
            network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
        isom_size += 1
    # TESTING FOR CORRECTNESS WHILE RUNNING
    #        if not Isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
    #            print("not unlabeled isom")
    #            print(seq_from_1)
    #            print(seq_from_2)
    #            print(isom_1_2)
    #            print(network1.subgraph(isom_1_2.keys()).edges())
    #            print(network2.subgraph(isom_2_1.keys()).edges())

    # Add the root to the isomorphism, if it was there
    isom_1_2[root1] = root2
    isom_2_1[root2] = root1

    # invert seq_from_2, rename to node names of network1, and append to seq_from_1
    return seq_from_1 + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)


# Code used for cases 1 and 2 in the head move red line algorithm
# returns moves1,moves2,addedNode1,addedNode2
def RL_Case1(N1, N2, x_1, isom_N1_N2, isom_N2_N1, randomNodes=False):
    """
    An implementation of Algorithm 7. Finds a sequence of head moves that makes it possible to add the highest tree node x_1 to the up-closed isomrophism.

    :param N1: a phylogenetic network.
    :param N2: a phylogenetic network.
    :param x_1: a tree node in N1.
    :param isom_N1_N2: a dictionary, containing a partial (up-closed) isomorphism map from N1 to N2. The inverse of isom_N2_N1.
    :param isom_N2_N1: a dictionary, containing a partial (up-closed) isomorphism map from N2 to N1. The inverse of isom_N1_N2.
    :param randomNodes: a boolean value, determining whether the random version of this algorithm is used.
    :return: a list of head moves in N1, a list of head moves in N2, a node of N1, a node of N2. After applying the moves to the networks, the nodes can be added to the up-closed isomorphism.
    """
    p_1 = Parent(N1, x_1)
    p_2 = isom_N1_N2[p_1]
    x_2 = Child(N2, p_2, exclude=isom_N2_N1.keys(), randomNodes=randomNodes)
    if N2.out_degree(x_2) == 2:
        #        print("Case tree")
        return [], [], x_1, x_2
    elif N2.in_degree(x_2) == 2:
        #        print("Case retic")
        c_2 = FindTreeNode(N2, excludedSet=isom_N2_N1.keys(), randomNodes=randomNodes)
        t_2 = Parent(N2, c_2)
        b_2 = Child(N2, c_2, exclude=[x_2], randomNodes=randomNodes)
        #        print(x_1,p_1,p_2,x_2,c_2,t_2,b_2)
        # TODO: Not in the latex algo, just a minor improvement:
        # If the other child $c_2$ of $p_2$ is a retic, then we can add it.
        if p_2 == t_2:
            return [], [], x_1, c_2
        if CheckMovable(N2, (p_2, x_2), x_2):
            q_2 = Parent(N2, x_2, exclude=[p_2])
            if x_2 == t_2:
                #                print("(p,x) movable and x=t")
                return [], [((q_2, x_2), x_2, (c_2, b_2))], x_1, c_2
            else:
                #                print("(p,x) movable and x!=t")
                return [], [((p_2, x_2), x_2, (t_2, c_2)), ((t_2, x_2), x_2, (c_2, b_2))], x_1, c_2
        else:
            d_2 = Child(N2, x_2)
            z_2 = Parent(N2, x_2, exclude=[p_2])
            # Find a leaf with parent not equal to d_2
            l_2 = FindLeaf(N2, excludedParents=[d_2], randomNodes=randomNodes)
            w_2 = Parent(N2, l_2)
            #            print("other")
            #            print(d_2,z_2,w_2,l_2)
            if l_2 == b_2:  # after the first move, b_2 is not a child of c_2 anymore, so we fix this by taking the other child of c_2 as b_2
                b_2 = Child(N2, c_2, exclude=[l_2, x_2], randomNodes=randomNodes)
            if d_2 == t_2:  # i.e., if, after the first move, x_2 is the parent of c_2
                return [], [((z_2, d_2), d_2, (w_2, l_2)), ((z_2, x_2), x_2, (c_2, b_2))], x_1, c_2
            else:
                return [], [((z_2, d_2), d_2, (w_2, l_2)), ((p_2, x_2), x_2, (t_2, c_2)),
                            ((t_2, x_2), x_2, (c_2, b_2))], x_1, c_2
    else:
        #        print("Case leaf")
        c_2 = FindTreeNode(N2, excludedSet=isom_N2_N1.keys(), randomNodes=randomNodes)
        t_2 = Parent(N2, c_2)

        if p_2 == t_2:
            #            print("no move")
            return [], [], x_1, c_2

            # Find a reticulation arc (s_2,r_2) that can be moved to (p_2,x_2)
        # No randomness required, because this arc will end up at its original position again.
        s_2 = None
        r_2 = None
        for node in N2.nodes():
            if N2.in_degree(node) == 2:
                for parent in N2.predecessors(node):
                    if parent != p_2 and CheckMovable(N2, (parent, node), node):
                        s_2 = parent
                        r_2 = node
            if s_2 != None:
                break
        q_2 = Parent(N2, r_2, exclude=[s_2])
        w_2 = Child(N2, r_2)
        #        print(p_2,x_2,t_2,c_2,s_2,q_2,r_2,w_2)
        if r_2 == p_2:
            #            print("r=p")
            if s_2 == t_2:
                return [], [((q_2, p_2), p_2, (t_2, c_2))], x_1, c_2
            elif q_2 == t_2:
                return [], [((s_2, p_2), p_2, (t_2, c_2))], x_1, c_2
            else:
                return [], [((s_2, p_2), p_2, (t_2, c_2)), ((t_2, p_2), p_2, (q_2, x_2)),
                            ((q_2, p_2), p_2, (s_2, c_2))], x_1, c_2
        if r_2 == t_2:
            #            print("r=t")
            if p_2 == q_2:
                return [], [((s_2, r_2), r_2, (p_2, x_2))], x_1, c_2
            else:
                return [], [((s_2, r_2), r_2, (p_2, x_2)), ((p_2, r_2), r_2, (q_2, c_2)),
                            ((q_2, r_2), r_2, (s_2, x_2))], x_1, c_2
        if s_2 != t_2:
            #            print("s!=t")
            return [], [((s_2, r_2), r_2, (p_2, x_2)), ((p_2, r_2), r_2, (t_2, c_2)), ((t_2, r_2), r_2, (s_2, x_2)),
                        ((s_2, r_2), r_2, (q_2, w_2))], x_1, c_2
        else:
            #            print("other")
            return [], [((s_2, r_2), r_2, (p_2, x_2)), ((p_2, r_2), r_2, (s_2, c_2)),
                        ((s_2, r_2), r_2, (q_2, w_2))], x_1, c_2

        # Code used for case 3 in the head move red line algorithm


# returns moves1,moves2,addedNode1,addedNode2
def RL_Case3(N1, N2, x_1, isom_N1_N2, isom_N2_N1, randomNodes=False):
    """
    An implementation of Algorithm 8. Finds a sequence of head moves that makes it possible to add the highest reticulation x_1 to the up-closed isomrophism.

    :param N1: a phylogenetic network.
    :param N2: a phylogenetic network.
    :param x_1: a highest reticulation node in N1.
    :param isom_N1_N2: a dictionary, containing a partial (up-closed) isomorphism map from N1 to N2. The inverse of isom_N2_N1.
    :param isom_N2_N1: a dictionary, containing a partial (up-closed) isomorphism map from N2 to N1. The inverse of isom_N1_N2.
    :param randomNodes: a boolean value, determining whether the random version of this algorithm is used.
    :return: a list of head moves in N1, a list of head moves in N2, a node of N1, a node of N2. After applying the moves to the networks, the nodes can be added to the up-closed isomorphism.
    """
    p_1 = Parent(N1, x_1)
    q_1 = Parent(N1, x_1, exclude=[p_1])
    p_2 = isom_N1_N2[p_1]
    cp_2 = Child(N2, p_2, exclude=isom_N2_N1.keys(), randomNodes=randomNodes)
    q_2 = isom_N1_N2[q_1]
    cq_2 = Child(N2, q_2, exclude=isom_N2_N1.keys(), randomNodes=randomNodes)

    # at least one tree node, so give up
    if N2.out_degree(cp_2) == 2 or N2.out_degree(cq_2) == 2:
        # The proof does not provide a sequence when at least of the nodes cp_2 or cq_2 is a tree node
        # TODO: If one of (p_2,cp_2) or (q_2,cq_2) is movable, we can still do something quite similar to what follows in the last case
        return [], [], None, None
    # Case 3ai
    if cp_2 == cq_2:
        #        print("no move")
        return [], [], x_1, cp_2
    # Case 3av
    elif N2.out_degree(cp_2) == 0 and N2.out_degree(cq_2) == 0:
        #        print("both leaves")
        # Find a head-movable arc (s_2,r_2)
        # As each reticulation node has a movable arc, we pick a random reticulation, and then a random movable arc incident to this node
        r_2 = FindRetic(N2, excludedSet=isom_N2_N1.keys(), randomNodes=randomNodes)
        s_2 = None
        for parent in N2.predecessors(r_2):
            if CheckMovable(N2, (parent, r_2), r_2):
                if not randomNodes:
                    s_2 = parent
                    break
                elif s_2 == None or random.getrandbits(1):
                    s_2 = parent
        if s_2 == p_2:
            #            print("s=p")
            return [], [((s_2, r_2), r_2, (q_2, cq_2)), ((q_2, r_2), r_2, (p_2, cp_2))], x_1, r_2
        else:
            #            print("s!=p")
            return [], [((s_2, r_2), r_2, (p_2, cp_2)), ((p_2, r_2), r_2, (q_2, cq_2))], x_1, r_2
    # Cases 3a(ii, iii, iv)
    else:
        #        print("at least one non-leaf")
        if nx.has_path(N2, cq_2, cp_2) or N2.out_degree(cp_2) == 0:
            # Swap p and q
            q_2, p_2 = p_2, q_2
            cp_2, cq_2 = cq_2, cp_2
        #        print(p_2 , q_2, cp_2, cq_2)
        if CheckMovable(N2, (p_2, cp_2), cp_2):
            #            print("movable cp_2")
            return [], [((p_2, cp_2), cp_2, (q_2, cq_2))], x_1, cp_2
        else:
            #            print("immovable cp_2")
            z = Child(N2, cp_2)
            t = Parent(N2, cp_2, exclude=[p_2])
            if t == q_2:
                return [], [], x_1, cp_2
            return [], [((cp_2, z), z, (q_2, cq_2)), ((t, cp_2), cp_2, (z, cq_2))], x_1, z


# Permutes the leaves of a network so that it becomes leaf isomorphic, provided the networks were already (non-labeled) ismorphic
def Permute_Leaves_Head(network1, network2, isom_1_2, isom_2_1, label_dict_1, label_dict_2):
    """
    An implementation of Algorithm 9. Determines a sequence of head moves that makes two isomorphic networks labeled isomorphic.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param isom_1_2: a dictionary, containing a partial (up-closed) isomorphism map from network1 to network2. The inverse of isom_2_1.
    :param isom_2_1: a dictionary, containing a partial (up-closed) isomorphism map from network2 to network1. The inverse of isom_1_2.
    :param label_dict_1: a dictionary, containing the correspondence of nodes of network1 with the leaf labels: keys are labels and values are nodes.
    :param label_dict_2: a dictionary, containing the correspondence of nodes of network2 with the leaf labels: keys are labels and values are nodes.
    :return: a list of head moves that turns network1 into network2.
    """
    #    for i,k in isom_1_2.items():
    #        print(i,k)
    sequence = []
    # labeldict[label]=leaf
    Y = list(label_dict_1.values())
    cycles = []
    while len(Y) > 0:
        y1_1 = Y.pop()
        y_2 = isom_1_2[y1_1]
        cycle = [y1_1]
        while network2.node[y_2].get('label') != network1.node[cycle[0]].get('label'):
            y_new_1 = label_dict_1[network2.node[y_2].get('label')]
            #            if len(set(network1.predecessors(cycle[-1]))&set(network1.predecessors(y_new_1)))==0:#cycle[-1] and y_new_1 have NO common parent
            #                cycle+=[y_new_1]
            cycle += [y_new_1]

            y_2 = isom_1_2[y_new_1]
            Y.remove(y_new_1)
        if len(cycle) > 1:
            cycles += [list(reversed(cycle))]
    #    print("cycles",cycles)

    t = None
    r = None
    for node in network1:
        if network1.in_degree(node) == 2:
            for parent in network1.predecessors(node):
                if CheckMovable(network1, (parent, node), node):
                    t = parent
                    r = node
        if r:
            break
    c_last = Child(network1, r)  # In the proof: c', the child of r in N
    s_last = Parent(network1, r, exclude=[t])

    for cycle in cycles:
        c = cycle

        # shift the cycle by one of t is the parent of the last leaf in the cycle
        if t in network1.predecessors(cycle[-1]):
            c = [cycle[-1]] + cycle[:-1]

        # Skip the first move if the head r of the moving arc is already above the last leaf in the cycle
        p_last = Parent(network1, c[-1])
        if r != p_last:
            move = ((t, r), r, (p_last, c[-1]))
            sequence.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
            network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
        moved_arc = (t, r)

        c_last_before = c_last

        for i in reversed(range(len(c))):
            p_i = Parent(network1, r, exclude=[moved_arc[0]])
            p_im1 = Parent(network1, c[i - 1])
            # if p_i==p_im1, then two consecutive leaves are in a cherry, and we can skip one of them.
            if p_i == p_im1:
                continue
            #                print("do nothing, swapping a cherry")
            else:
                move = ((p_i, r), r, (p_im1, c[i - 1]))
                sequence.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
                network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
                moved_arc = (p_i, r)
            # If c in the `original position' (s,c) of (t,r) is permuted by the cycle, then we need to change this original position accordingly
            if c[i] == c_last_before:
                c_last = c[i - 1]

    # Put the moving arc back to its original position
    # The proof does this for every cycle, but can also be done once at the end
    if (s_last, c_last) in network1.edges():
        # (t,r) might already be back at this position
        move = ((t, r), r, (s_last, c_last))
        sequence.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
        network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)

    #    print("isomorphic", Isomorphic(network1,network2))
    #    for s in sequence:
    #        print(s[0])
    return sequence


# A heuristic for finding a head move sequence between two networks.
def Red_Line(network1, network2):
    """
    An implementation of Algorithm 10. Finds a sequence of head moves from network1 to network2 by building an up-closed isomorphism.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :return: A list of head moves from network1 to network2.
    """
    # Find the root and labels of the networks
    root1 = Root(network1)
    root2 = Root(network2)
    label_dict_1 = Labels(network1)
    label_dict_2 = Labels(network2)

    # initialize isomorphism
    isom_1_2 = dict()
    isom_1_2[root1] = root2
    isom_2_1 = dict()
    isom_2_1[root2] = root1
    isom_size = 1

    # Check if the roots are of the same type
    if network1.out_degree(root1) != network2.out_degree(root2):
        return False

    # Keep track of the size of the isomorphism and the size it is at the end of the red line algorithm
    goal_size = len(network1) - len(label_dict_1)

    # init lists of sequence of moves
    # list of (moving_edge,moving_endpoint,from_edge,to_edge)
    seq_from_1 = []
    seq_from_2 = []
    # TODO keep track of highest nodes?

    # Do the red line algorithm
    while (isom_size < goal_size):
        highest_tree_node_network1, highest_retic_network1 = HighestNodesBelow(network1, isom_1_2.keys())
        highest_tree_node_network2, highest_retic_network2 = HighestNodesBelow(network2, isom_2_1.keys())

        # Case1
        if highest_tree_node_network1 != None:
            moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case1(network1, network2,
                                                                                                  highest_tree_node_network1,
                                                                                                  isom_1_2, isom_2_1)
        # Case2
        elif highest_tree_node_network2 != None:
            moves_network_2, moves_network_1, added_node_network2, added_node_network1 = RL_Case1(network2, network1,
                                                                                                  highest_tree_node_network2,
                                                                                                  isom_2_1, isom_1_2)
        # Case3
        else:
            moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case3(network1, network2,
                                                                                                  highest_retic_network1,
                                                                                                  isom_1_2, isom_2_1)

        # Now perform the moves and update the isomorphism
        isom_1_2[added_node_network1] = added_node_network2
        isom_2_1[added_node_network2] = added_node_network1
        for move in moves_network_1:
            seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
            network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
        for move in moves_network_2:
            seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
            network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
        isom_size += 1
    # TESTING FOR CORRECTNESS WHILE RUNNING
    #        if not nx.is_isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
    #            print("not unlabeled isom")
    #            print(seq_from_1)
    #            print(seq_from_2)
    #            print(network1.subgraph(isom_1_2.keys()).edges())
    #            print(network2.subgraph(isom_2_1.keys()).edges())

    # TODO Debugging, remove after for speed
    #    if not nx.is_isomorphic(network1,network2):
    #        print("not unlabeled isom")
    #        print(network1.edges())
    #        print(network2.edges())
    #    else:
    #        print("unlabeled isomorphic :)")

    # Add the leaves to the isomorphism
    for node_1 in network1.nodes():
        if network1.out_degree(node_1) == 0:
            parent_1 = Parent(network1, node_1)
            parent_2 = isom_1_2[parent_1]
            node_2 = Child(network2, parent_2, exclude=isom_2_1.keys())
            isom_1_2[node_1] = node_2
            isom_2_1[node_2] = node_1

    # Permute the leaves
    seq_permute = Permute_Leaves_Head(network1, network2, isom_1_2, isom_2_1, label_dict_1, label_dict_2)

    # invert seq_from_2, rename to node names of network1, and append to seq_from_1
    return seq_from_1 + seq_permute + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)


def Red_Line_Random(network1, network2, repeats=1):
    """
    Finds a sequence of head moves from network1 to network2 by randomly building an up-closed isomorphism a number of times, and only keeping the shortest sequence.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :param repeats: an integer, determining how many repeats of Red_Line_Random_Single are performed.
    :return: A list of head moves from network1 to network2.
    """
    best_seq = None
    for i in range(repeats):
        candidate_seq = Red_Line_Random_Single(network1, network2)
        if best_seq == None or len(best_seq) > len(candidate_seq):
            best_seq = candidate_seq
    return best_seq


def Red_Line_Random_Single(network1, network2):
    """
    An implementation of Algorithm 11. Finds a sequence of head moves from network1 to network2 by building randomly an up-closed isomorphism.

    :param network1: a phylogenetic network.
    :param network2: a phylogenetic network.
    :return: A list of head moves from network1 to network2.
    """
    # Find the root and labels of the networks
    root1 = Root(network1)
    root2 = Root(network2)
    label_dict_1 = Labels(network1)
    label_dict_2 = Labels(network2)

    # initialize isomorphism
    isom_1_2 = dict()
    isom_1_2[root1] = root2
    isom_2_1 = dict()
    isom_2_1[root2] = root1
    isom_size = 1

    # Check if the roots are of the same type
    if network1.out_degree(root1) != network2.out_degree(root2):
        return False

    # Keep track of the size of the isomorphism and the size it is at the end of the red line algorithm
    goal_size = len(network1) - len(label_dict_1)

    # init lists of sequence of moves
    # list of (moving_edge,moving_endpoint,from_edge,to_edge)
    seq_from_1 = []
    seq_from_2 = []
    # TODO keep track of highest nodes?

    # Do the red line algorithm
    while (isom_size < goal_size):
        highest_tree_node_network1, highest_retic_network1 = HighestNodesBelow(network1, isom_1_2.keys(), allnodes=True)
        highest_tree_node_network2, highest_retic_network2 = HighestNodesBelow(network2, isom_2_1.keys(), allnodes=True)

        # Construct a list of all highest nodes in a tuple with the corresponding network (in random order)
        # I.e. If u is a highest node of network one, it will appear in the list as (u,1)
        highest_nodes_network1 = [(u, 1) for u in highest_tree_node_network1 + highest_retic_network1]
        highest_nodes_network2 = [(u, 2) for u in highest_tree_node_network2 + highest_retic_network2]
        candidate_highest_nodes = highest_nodes_network1 + highest_nodes_network2
        random.shuffle(candidate_highest_nodes)

        # As some cases do not give an addition to the isom, we continue trying lowest nodes until we find one that does.
        for highest_node, network_number in candidate_highest_nodes:
            u = highest_node
            # Case1
            if network_number == 1 and network1.in_degree(u) == 1:
                #                print("Case 1")
                moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case1(network1,
                                                                                                      network2, u,
                                                                                                      isom_1_2,
                                                                                                      isom_2_1,
                                                                                                      randomNodes=True)
                break
            # Case2
            elif network_number == 2 and network2.in_degree(u) == 1:
                #                print("Case 2")
                moves_network_2, moves_network_1, added_node_network2, added_node_network1 = RL_Case1(network2,
                                                                                                      network1, u,
                                                                                                      isom_2_1,
                                                                                                      isom_1_2,
                                                                                                      randomNodes=True)
                break
            # Case3
            elif network_number == 1 and network1.in_degree(u) == 2:
                #                print("Case 3")
                moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case3(network1,
                                                                                                      network2, u,
                                                                                                      isom_1_2,
                                                                                                      isom_2_1,
                                                                                                      randomNodes=True)
                if added_node_network2:
                    break
            # Case3'
            elif network_number == 2 and network2.in_degree(u) == 2:
                #                print("Case 3'")
                moves_network_2, moves_network_1, added_node_network2, added_node_network1 = RL_Case3(network2,
                                                                                                      network1, u,
                                                                                                      isom_2_1,
                                                                                                      isom_1_2,
                                                                                                      randomNodes=True)
                if added_node_network2:
                    break

        # Now perform the moves and update the isomorphism
        isom_1_2[added_node_network1] = added_node_network2
        isom_2_1[added_node_network2] = added_node_network1
        for move in moves_network_1:
            seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
            network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
        for move in moves_network_2:
            seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
            network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
        isom_size += 1
    # TESTING FOR CORRECTNESS WHILE RUNNING
    #        if not nx.is_isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
    #            print("not unlabeled isom")
    #            print(seq_from_1)
    #            print(seq_from_2)
    #            print(network1.subgraph(isom_1_2.keys()).edges())
    #            print(network2.subgraph(isom_2_1.keys()).edges())
    #            print(network1.edges())
    #            print(network2.edges())

    # TODO Debugging, remove after for speed
    #    if not nx.is_isomorphic(network1,network2):
    #        print("not unlabeled isom")
    #        print(network1.edges())
    #        print(network2.edges())
    #    else:
    #        print("unlabeled isomorphic :)")

    # Add the leaves to the isomorphism
    for node_1 in network1.nodes():
        if network1.out_degree(node_1) == 0:
            parent_1 = Parent(network1, node_1)
            parent_2 = isom_1_2[parent_1]
            node_2 = Child(network2, parent_2, exclude=isom_2_1.keys())
            isom_1_2[node_1] = node_2
            isom_2_1[node_2] = node_1

    # Permute the leaves
    seq_permute = Permute_Leaves_Head(network1, network2, isom_1_2, isom_2_1, label_dict_1, label_dict_2)

    # invert seq_from_2, rename to node names of network1, and append to seq_from_1
    return seq_from_1 + seq_permute + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)


# Find the original location of the moving_endpoint,
# That is, the edge from which we remove it.
def From_Edge(network, moving_edge, moving_endpoint):
    """
    Finds the original location (from-edge) of the moving_endpoint of an edge if it is moved.

    :param network: a phylogenetic network.
    :param moving_edge: an edge of the network.
    :param moving_endpoint: a node of the network, which must be an endpoint of the edge.
    :return: a pair of nodes (p,c) where p and c are a parent and child of the moving_endpoint such that they are both not part of the moving_edge.
    """
    other_parent = Parent(network, moving_endpoint, exclude=moving_edge)
    other_child = Child(network, moving_endpoint, exclude=moving_edge)
    return (other_parent, other_child)


def InvertMoveSequence(seq):
    """
    Inverts a sequence of moves.

    :param seq: a list of moves, where each moves is in the format (moving_edge,moving_endpoint,from_edge,to_edge).
    :return: the inverse list of moves.
    """
    newSeq = []
    for move in reversed(seq):
        moving_edge, moving_endpoint, from_edge, to_edge = move
        newSeq.append((moving_edge, moving_endpoint, to_edge, from_edge))
    return newSeq


def ReplaceNodeNamesInMoveSequence(seq, isom):
    """
    Renames the nodes in a sequence of moves using an isomorphism mapping between two networks.

    :param seq: a list of moves, implicitly using the nodes of a network.
    :param isom: a dictionary, containing a bijective mapping from nodes of the networks to another set.
    :return: a list of moves where the nodes are replaced by their image under the isom mapping.
    """
    if type(seq) == int:
        return isom[seq]
    return list(map(lambda x: ReplaceNodeNamesInMoveSequence(x, isom), seq))


def ReplaceNodeNamesByOriginal(network, seq):
    """
    Renames the nodes in a sequence of moves by their original names as given in the input. These are stored as a node attribute 'original'.

    :param network: a phylogenetic network.
    :param seq: a sequence of moves on teh network.
    :return: a list of moves on the network using the original node names as used in the input.
    """
    if type(seq) == int:
        return network.node[seq]['original']
    if seq == 'rho':
        return "rho"
    return list(map(lambda x: ReplaceNodeNamesByOriginal(network, x), seq))


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
