# # Returns all nodes below a given node (including the node itself)
# def AllBelow(network, node):
#     """
#     Finds all nodes below a given node in a network.

#     :param network: a phylogenetic network.
#     :param node: a node in the network.
#     :return: all nodes in the network that are below the chosen node, including this node.
#     """
#     lengths = nx.single_source_shortest_path_length(network, node)
#     return lengths.keys()


# # Find the lowest nodes above a given set.
# # The excluded set must include all leaves.
# def LowestReticAndTreeNodeAbove(network, excludedSet, allnodes=False):
#     """
#     Finds a list of lowest tree nodes and a list of lowest reticulation nodes above a given set of nodes.

#     :param network: a phylogenetic network.
#     :param excludedSet: a set of nodes of the network.
#     :param allnodes: a boolean value that determines whether we try to find all lowest nodes (True) or only one lowest node of each type (False, Default).
#     :return: A list of tree nodes and a list of reticulation nodes, so that each element of these lists has all their children in the excludedSet. If not allNodes, then these lists have length at most 1.
#     """
#     retic = None
#     tree_node = None
#     lowest_retics = []
#     lowest_tree_nodes = []
#     for node in network.nodes():
#         if node not in excludedSet:
#             for c in network.successors(node):
#                 if c not in excludedSet:
#                     break
#             # else runs if the loop was not ended by a break
#             # this happens exactly when all of the children are in excludedSet
#             else:
#                 if network.out_degree(node) == 2:
#                     if allnodes:
#                         lowest_tree_nodes += [node]
#                     elif tree_node == None:
#                         # For simplicity in description, take the FIRST lowest tree node that we encounter (sim. for the reticulations)
#                         tree_node = node
#                 elif network.in_degree(node) == 2:
#                     if allnodes:
#                         lowest_retics += [node]
#                     elif retic == None:
#                         retic = node
#                 if not allnodes and tree_node != None and retic != None:
#                     # stop if both types of lowest nodes are found
#                     break
#     if allnodes:
#         return lowest_tree_nodes, lowest_retics
#     return tree_node, retic


# # Find the highest nodes below a given set.
# # The excluded set must include the root.
# def HighestNodesBelow(network, excludedSet, allnodes=False):
#     """
#     Finds a list of highest tree nodes and a list of highest reticulation nodes below a given set of nodes.

#     :param network: a phylogenetic network.
#     :param excludedSet: a set of nodes of the network.
#     :param allnodes: a boolean value that determines whether we try to find all highest nodes (True) or only one highest node of each type (False, Default).
#     :return: A list of tree nodes and a list of reticulation nodes, so that each element of these lists has all their parents in the excludedSet. If not allNodes, then these lists have length at most 1.
#     """
#     retic = None
#     tree_node = None
#     highest_retics = []
#     highest_tree_nodes = []
#     for node in network.nodes():
#         if node not in excludedSet:
#             for c in network.predecessors(node):
#                 if c not in excludedSet:
#                     break
#             # else runs if the loop was not ended by a break
#             # this happens exactly when all of the parents are in excludedSet
#             else:
#                 if network.out_degree(node) == 2:
#                     if allnodes:
#                         highest_tree_nodes += [node]
#                     elif tree_node == None:
#                         # For simplicity in description, take the FIRST highest tree node that we encounter (sim. for the reticulations and leaves)
#                         tree_node = node
#                 elif network.in_degree(node) == 2:
#                     if allnodes:
#                         highest_retics += [node]
#                     elif retic == None:
#                         retic = node
#                 if not allnodes and retic != None and tree_node != None:
#                     # stop if all types of highest nodes are found
#                     break
#     if allnodes:
#         return highest_tree_nodes, highest_retics
#     return tree_node, retic


# def FindTreeNode(network, excludedSet=[], randomNodes=False):
#     """
#     Finds a (random) tree node in a network.

#     :param network: a phylogenetic network.
#     :param excludedSet: a set of nodes of the network.
#     :param randomNodes: a boolean value.
#     :return: a tree node of the network not in the excludedSet, or None if no such node exists. If randomNodes, then a tree node is selected from all candidates uniformly at random.
#     """
#     all_found = []
#     for node in network.nodes():
#         if node not in excludedSet and network.out_degree(node) == 2 and network.in_degree(node) == 1:
#             if randomNodes:
#                 all_found += [node]
#             else:
#                 return node
#     if all_found and randomNodes:
#         return random.choice(all_found)
#     return None


# def FindLeaf(network, excludedSet=[], excludedParents=[], randomNodes=False):
#     """
#     Finds a (random) leaf in a network.

#     :param network: a phylogenetic network.
#     :param excludedSet: a set of nodes of the network.
#     :param excludedParents: a set of nodes of the network.
#     :param randomNodes: a boolean value.
#     :return: a leaf of the network not in the excludedSet so that its parent is nt in excludedParents, or None if no such node exists. If randomNodes, then a leaf is selected from all candidates uniformly at random.
#     """
#     all_found = []
#     for node in network.nodes():
#         parent = Parent(network, node)
#         if network.out_degree(node) == 0 and parent not in excludedParents and node not in excludedSet:
#             if randomNodes:
#                 all_found += [node]
#             else:
#                 return node
#     if all_found and randomNodes:
#         return random.choice(all_found)
#     return None


# def FindRetic(network, excludedSet=[], randomNodes=False):
#     """
#     Finds a (random) reticulation in a network.

#     :param network: a phylogenetic network.
#     :param excludedSet: a set of nodes of the network.
#     :param randomNodes: a boolean value.
#     :return: a reticulation node of the network not in the excludedSet, or None if no such node exists. If randomNodes, then a reticulation is selected from all candidates uniformly at random.
#     """
#     all_found = []
#     for node in network.nodes():
#         if node not in excludedSet and network.in_degree(node) == 2:
#             if randomNodes:
#                 all_found += [node]
#             else:
#                 return node
#     if all_found and randomNodes:
#         return random.choice(all_found)
#     return None


# def Parent(network, node, exclude=[], randomNodes=False):
#     """
#     Finds a parent of a node in a network.

#     :param network: a phylogenetic network.
#     :param node: a node in the network.
#     :param exclude: a set of nodes of the network.
#     :param randomNodes: a boolean value.
#     :return: a parent of node that is not in the set of nodes exclude. If randomNodes, then this parent is selected uniformly at random from all candidates.
#     """
#     parent = None
#     for p in network.predecessors(node):
#         if p not in exclude:
#             if not randomNodes:
#                 return p
#             elif parent == None or random.getrandbits(1):
#                 # As there are at most two parents, we can simply replace the previous parent with probability .5 to get a random parent
#                 parent = p
#     return parent


# def Child(network, node, exclude=[], randomNodes=False):
#     """
#     Finds a child node of a node in a network.

#     :param network: a phylogenetic network.
#     :param node: a node in the network.
#     :param exclude: a set of nodes of the network.
#     :param randomNodes: a boolean value.
#     :return: a child of node that is not in the set of nodes exclude. If randomNodes, then this child node is selected uniformly at random from all candidates.
#     """
#     child = None
#     for c in network.successors(node):
#         if c not in exclude:
#             if not randomNodes:
#                 return c
#             elif child == None or random.getrandbits(1):
#                 # As there are at most two children, we can simply replace the previous child with probability .5 to get a random parent
#                 child = c
#     return child


# # Returns the root of a network
# def Root(network):
#     """
#     Finds the root of a phylogenetic network.

#     :param network: a phylogenetic network.
#     :return: the root node of this network.
#     """
#     for node in network.nodes():
#         if network.in_degree(node) == 0:
#             return node
#     return None


# # Returns a dictionary with node labels, keyed by the labels
# def Labels(network):
#     """
#     Returns the correspondence between the leaves and the leaf-labels of a given network

#     :param network: a phylogenetic network
#     :return: a dictionary, where the keys are labels and the values are nodes of the network.
#     """
#     label_dict = dict()
#     for node in network.nodes():
#         node_label = network.node[node].get('label')
#         if node_label:
#             label_dict[node_label] = node
#     return label_dict


#     # Find a sequence by choosing the move that most decreases the upper bound on the number of moves
# # This works as long as we can always decrease the bound.
# # E.g.1, this upper bound can be the length of the sequence given by Green_Line(N1,N2), the bound can always decrease after one move, if we take the move from the GL sequence (IMPLEMENTED)
# # TODO E.g.2, take the upper bound given by this algorithm with bound Green_Line
# def Deep_Dive_Scored(network1, network2, head_moves=True, bound_heuristic=Green_Line):
#     """
#     An experimental method that returns a sequence of tail/rSPR moves from network1 to network2, using the isomorphism-building heuristic for the chosen type of moves.


#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :param head_moves: a boolean value that determines whether head moves are used in addition to tail moves. If True we use rSPR moves, if False we use only tail moves.
#     :param bound_heuristic: a heuristic that finds a sequence between the networks quickly.
#     :return: a sequence of moves from network1 to network2.
#     """
#     if Isomorphic(network1, network2):
#         return []
#     seq = []
#     current_network = network1
#     current_best = []
#     for move in bound_heuristic(network1, network2, head_moves=head_moves):
#         current_best += [(move[0], move[1], move[3])]
#     if not current_best:
#         return False
#     done = False
#     current_length = 0
#     while not done:
#         candidate_moves = AllValidMoves(current_network, tail_moves=True, head_moves=head_moves)
#         for move in candidate_moves:
#             candidate_network = DoMove(current_network, *move)
#             if Isomorphic(candidate_network, network2):
#                 return current_best[:current_length] + [move]
#             candidate_sequence = bound_heuristic(candidate_network, network2, head_moves=head_moves)
#             if current_length + len(candidate_sequence) + 1 < len(current_best):
#                 current_best = current_best[:current_length] + [move]
#                 for move2 in candidate_sequence:
#                     current_best += [(move2[0], move2[1], move2[3])]
#         next_move = current_best[current_length]
#         current_network = DoMove(current_network, *next_move)
#         current_length += 1
#     return True
