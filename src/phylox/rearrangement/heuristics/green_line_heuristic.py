# def GL_Case1_rSPR(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False):
#     """
#     An implementation of Algorithm 2. Finds a sequence of rSPR moves that makes it possible to add the lowest reticulation up to the down-closed isomrophism.

#     :param N: a phylogenetic network.
#     :param Np: a phylogenetic network.
#     :param up: a lowest reticulation node of Np above the isomorphism.
#     :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
#     :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
#     :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
#     :return: a list of rSPR moves in N, a list of rSPR moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism.
#     """
#     # use notation as in the paper
#     # ' is denoted p
#     xp = Child(Np, up)
#     x = isom_Np_N[xp]
#     z = Parent(N, x, exclude=isom_N_Np.keys(), randomNodes=randomNodes)

#     # Case1a: z is a reticulation
#     if N.in_degree(z) == 2:
#         return [], [], z, up
#     # Case1b: z is not a reticulation
#     # Find a retic v in N not in the isom yet
#     v = FindRetic(N, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes)
#     u = None
#     for parent in N.predecessors(v):
#         if CheckValid(N, (parent, v), v, (z, x)):
#             if not randomNodes:
#                 u = parent
#                 break
#             elif u == None or random.getrandbits(1):
#                 u = parent
#     # If v has an incoming arc (u,v) movable to (z,x), perform this move
#     if u != None:
#         return [((u, v), v, (z, x))], [], v, up
#     # if none of the moves is valid, this must be because
#     # v is already a reticulation above x with its movable incoming arc (u,v) with u=z.
#     return [], [], v, up


# def GL_Case1_Tail(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False):
#     """
#     An implementation of Algorithm 6. Finds a sequence of tail moves that makes it possible to add the lowest reticulation up to the down-closed isomrophism.

#     :param N: a phylogenetic network.
#     :param Np: a phylogenetic network.
#     :param up: a lowest reticulation node of Np above the isomorphism.
#     :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
#     :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
#     :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
#     :return: a list of tail moves in N, a list of tail moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism. Returns false if the networks are not isomorphic with 2 leaves and 1 reticulation.
#     """
#     # use notation as in the paper
#     # ' is denoted p
#     xp = Child(Np, up)
#     x = isom_Np_N[xp]
#     z = Parent(N, x, exclude=isom_N_Np.keys(), randomNodes=randomNodes)
#     # Case1a: z is a reticulation
#     if N.in_degree(z) == 2:
#         return [], [], z, up
#     # Case1b: z is not a reticulation
#     # z is a tree node
#     if check_movable(N, (z, x), z):
#         # Case1bi: (z,x) is movable
#         # Find a reticulation u in N not in the isomorphism yet
#         # TODO: Can first check if the other parent of x suffices here, should heuristcally be better
#         u = FindRetic(N, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes)
#         v = Child(N, u)
#         if v == x:
#             return [], [], u, up
#         # we may now assume v!=x
#         if z == v:
#             v = Child(N, z, exclude=[x])
#             w = Parent(N, u, randomNodes=randomNodes)
#             return [((z, v), z, (w, u))], [], u, up
#         w = Parent(N, u, exclude=[z], randomNodes=randomNodes)
#         return [((z, x), z, (u, v)), ((z, v), z, (w, u))], [], u, up
#     # Case1bii: (z,x) is not movable
#     c = Parent(N, z)
#     d = Child(N, z, exclude=[x])
#     # TODO: b does not have to exist if we have an outdeg-2 root, this could be c!
#     b = Parent(N, c)
#     if N.in_degree(b) != 0:
#         # Case1biiA: b is not the root of N
#         a = Parent(N, b, randomNodes=randomNodes)
#         # First do the move ((c,d),c,(a,b)), then Case1bi applies as (z,x) is now movable
#         newN = DoMove(N, (c, d), c, (a, b))
#         u = FindRetic(newN, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes)
#         v = Child(newN, u)
#         if v == x:
#             # In this case, u is a reticulation parent of x and u is not in the isom. Hence, we can simply add it to the isom.
#             # Note: The tail move we did is not necessary!
#             # TODO: First check both parents of x for being a reticulation not in the isomorphism yet
#             return [], [], u, up
#         # we may now assume v!=x
#         if z == v:
#             # This only happens if z==v and u==b
#             # we move z up above the retic b as well, too
#             w = Parent(newN, b, randomNodes=randomNodes)
#             return [((c, d), c, (a, b)), ((z, d), z, (w, b))], [], u, up
#         w = Parent(newN, u, exclude=[z], randomNodes=randomNodes)
#         return [((c, d), c, (a, b)), ((z, x), z, (u, v)), ((z, v), z, (w, u))], [], u, up
#     # Case1biiB: b is the root of N
#     # Note: d is not in the isomorphism
#     e = Child(N, d)
#     if e == x:
#         return [], [], d, up
#     if N.out_degree(x) == 2:
#         s = Child(N, x)
#         t = Child(N, x, exclude=[s])
#         if s == e:
#             return [((x, t), x, (d, e))], [], d, up
#         if t == e:
#             return [((x, s), x, (d, e))], [], d, up
#         return [((x, s), x, (d, e)), ((x, e), x, (z, t)), ((x, t), x, (d, s))], [], d, up
#     if N.out_degree(e) == 2:
#         s = Child(N, e)
#         t = Child(N, e, exclude=[s])
#         if s == x:
#             return [((e, t), e, (z, x))], [], d, up
#         if t == x:
#             return [((e, s), e, (z, x))], [], d, up
#         return [((e, s), e, (z, x)), ((e, x), e, (d, t)), ((e, t), e, (z, s))], [], d, up
#     # neither are tree nodes, so both must be leaves
#     # In that case, there is no sequence between the two networks.
#     return [], [], None, None


# def GL_Case3(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False):
#     """
#     An implementation of Algorithm 3. Finds a sequence of tail moves that makes it possible to add the lowest tree node up to the down-closed isomrophism.

#     :param N: a phylogenetic network.
#     :param Np: a phylogenetic network.
#     :param up: a lowest tree node of Np above the isomorphism.
#     :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
#     :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
#     :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
#     :return: a list of tail moves in N, a list of tail moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism.
#     """
#     # Find the children x' and y' of u'
#     xp, yp = list(Np.successors(up))
#     # Find the corresponding nodes x and y in N
#     x = isom_Np_N[xp]
#     y = isom_Np_N[yp]
#     # Find the set of common parents of x and y
#     parents_x = set(N.predecessors(x))
#     parents_y = set(N.predecessors(y))
#     common_parents = parents_x & parents_y
#     # Case3a: x and y have a common parent not in the isom
#     common_parents_not_Y = []
#     for parent in common_parents:
#         if parent not in isom_N_Np.keys():
#             # then parent can be mapped to up
#             common_parents_not_Y += [parent]
#             if not randomNodes:
#                 return [], [], parent, up
#     if common_parents_not_Y:
#         return [], [], random.choice(common_parents_not_Y), up

#     # Case3b: x and y do not have a common parent in the isomorphism
#     # For both, find a parent not in the isomorphism yet
#     # TODO: preferably make them tree nodes
#     z_x = Parent(N, x, exclude=isom_N_Np.keys(), randomNodes=randomNodes)
#     z_y = Parent(N, y, exclude=isom_N_Np.keys(), randomNodes=randomNodes)

#     # Case3bi: (z_x,x) is movable
#     if CheckValid(N, (z_x, x), z_x, (z_y, y)):
#         return [((z_x, x), z_x, (z_y, y))], [], z_x, up
#     # Case3bii: (z_y,y) is movable
#     if CheckValid(N, (z_y, y), z_y, (z_x, x)):
#         return [((z_y, y), z_y, (z_x, x))], [], z_y, up
#     # Case3biii: Neither (z_x,x) nor (z_y,y) is movable

#     if N.in_degree(z_x) == 2 or N.in_degree(z_y) == 2:
#         return [], [], None, None
#     # Both these parents are tree nodes.
#     # This happens always in the non-random version, as otherwise there is a lowest reticulation node and we would be in Case 1 or 2.

#     # As both nodes are tree nodes and the arcs immovable, both arcs hang of the side of a triangle.
#     # Find the top node of the triangle for z_x
#     c_x = Parent(N, z_x)
#     b_x = Parent(N, c_x)

#     # Find the top node of the triangle for z_y
#     c_y = Parent(N, z_y)
#     #    print()
#     #    print(z_y,N.edges())

#     b_y = Parent(N, c_y)

#     if N.in_degree(b_x) == 0:
#         # c_x is the child of the root
#         # c_x!=c_y, so c_y is not the child of the root
#         # swap the roles of x and y
#         x, y = y, x
#         z_x, z_y = z_y, z_x
#         b_x, b_y = b_y, b_x
#         c_x, c_y = c_y, c_x
#     # c_x is not the child of the root
#     # find a parent of b_x, and the bottom node of the triangle d_x
#     a_x = Parent(N, b_x, randomNodes=randomNodes)
#     d_x = Child(N, c_x, exclude=[z_x])
#     return [((c_x, d_x), c_x, (a_x, b_x)), ((z_x, x), z_x, (z_y, y))], [], z_x, up


# def Green_Line(network1, network2, head_moves=True):
#     """
#     An implementation of Algorithm 4 and its tail move counterpart. Finds a sequence of tail/rSPR moves from network1 to network2 by building a down-closed isomorphism.

#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :param head_moves: a boolean value determining whether head moves are allowed (if True we use rSPR moves, if False we only use tail moves).
#     :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
#     """
#     # Find the root and labels of the networks
#     root1 = Root(network1)
#     root2 = Root(network2)
#     label_dict_1 = Labels(network1)
#     label_dict_2 = Labels(network2)

#     # initialize isomorphism
#     isom_1_2 = dict()
#     isom_2_1 = dict()
#     isom_size = 0
#     for label, node1 in label_dict_1.items():
#         node2 = label_dict_2[label]
#         isom_1_2[node1] = node2
#         isom_2_1[node2] = node1
#         isom_size += 1

#     # Keep track of the size of the isomorphism and the size it is at the end of the green line algorithm
#     goal_size = len(network1) - 1

#     # init lists of sequence of moves
#     # list of (moving_edge,moving_endpoint,from_edge,to_edge)
#     seq_from_1 = []
#     seq_from_2 = []
#     # TODO keep track of lowest nodes? (Currently, the code does not do this, but it could speed up the code)

#     # Do the green line algorithm
#     while (isom_size < goal_size):
#         # Find lowest nodes above the isom in the networks:
#         lowest_tree_node_network1, lowest_retic_network1 = LowestReticAndTreeNodeAbove(network1, isom_1_2.keys())
#         lowest_tree_node_network2, lowest_retic_network2 = LowestReticAndTreeNodeAbove(network2, isom_2_1.keys())

#         ######################################
#         # Case1: a lowest retic in network1
#         if lowest_retic_network1 != None:
#             # use notation as in the paper network1 = N', network2 = N, where ' is denoted p
#             up = lowest_retic_network1
#             if head_moves:
#                 moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_rSPR(network2,
#                                                                                                              network1,
#                                                                                                              up,
#                                                                                                              isom_2_1,
#                                                                                                              isom_1_2)
#             else:
#                 moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_Tail(network2,
#                                                                                                              network1,
#                                                                                                              up,
#                                                                                                              isom_2_1,
#                                                                                                              isom_1_2)
#                 if added_node_network_1 == None:
#                     return False
#         ######################################
#         # Case2: a lowest retic in network2
#         elif lowest_retic_network2 != None:
#             # use notation as in the paper network2 = N', network1 = N, where ' is denoted p
#             up = lowest_retic_network2
#             if head_moves:
#                 moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_rSPR(network1,
#                                                                                                              network2,
#                                                                                                              up,
#                                                                                                              isom_1_2,
#                                                                                                              isom_2_1)
#             else:
#                 moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_Tail(network1,
#                                                                                                              network2,
#                                                                                                              up,
#                                                                                                              isom_1_2,
#                                                                                                              isom_2_1)
#                 if added_node_network_1 == None:
#                     return False

#                     ######################################
#         # Case3: a lowest tree node in network1
#         else:
#             # use notation as in the paper network1 = N, network2 = N'
#             up = lowest_tree_node_network2
#             moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case3(network1, network2,
#                                                                                                     up, isom_1_2,
#                                                                                                     isom_2_1)
#         # Now perform the moves and update the isomorphism
#         isom_1_2[added_node_network_1] = added_node_network_2
#         isom_2_1[added_node_network_2] = added_node_network_1
#         for move in moves_network_1:
#             seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
#             network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
#         for move in moves_network_2:
#             seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
#             network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
#         isom_size += 1
#     # TESTING FOR CORRECTNESS WHILE RUNNING
#     #        if not Isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
#     #            print("not unlabeled isom")
#     #            print(seq_from_1)
#     #            print(seq_from_2)
#     #            print(network1.subgraph(isom_1_2.keys()).edges())
#     #            print(network2.subgraph(isom_2_1.keys()).edges())

#     # Add the root to the isomorphism, if it was there
#     isom_1_2[root1] = root2
#     isom_2_1[root2] = root1
#     # invert seq_from_2, rename to node names of network1, and append to seq_from_1
#     return seq_from_1 + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)


# def Green_Line_Random(network1, network2, head_moves=True, repeats=1):
#     """
#     Finds a sequence of tail/rSPR moves from network1 to network2 by randomly building a down-closed isomorphism a number of times, and only keeping the shortest sequence.

#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :param head_moves: a boolean value determining whether head moves are allowed (if True we use rSPR moves, if False we only use tail moves).
#     :param repeats: an integer, determining how many repeats of Green_Line_Random_Single are performed.
#     :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
#     """
#     best_seq = None
#     for i in range(repeats):
#         candidate_seq = Green_Line_Random_Single(network1, network2, head_moves=head_moves)
#         if candidate_seq == False:
#             return False
#         if best_seq == None or len(best_seq) > len(candidate_seq):
#             best_seq = candidate_seq
#     return best_seq


# def Green_Line_Random_Single(network1, network2, head_moves=True):
#     """
#     An implementation of Algorithm 5 and its tail move counterpart. Finds a sequence of tail/rSPR moves from network1 to network2 by randomly building a down-closed isomorphism.

#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :param head_moves: a boolean value determining whether head moves are allowed (if True we use rSPR moves, if False we only use tail moves).
#     :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
#     """
#     # Find the root and labels of the networks
#     root1 = Root(network1)
#     root2 = Root(network2)
#     label_dict_1 = Labels(network1)
#     label_dict_2 = Labels(network2)

#     # initialize isomorphism
#     isom_1_2 = dict()
#     isom_2_1 = dict()
#     isom_size = 0
#     for label, node1 in label_dict_1.items():
#         node2 = label_dict_2[label]
#         isom_1_2[node1] = node2
#         isom_2_1[node2] = node1
#         isom_size += 1

#     # Keep track of the size of the isomorphism and the size it is at the end of the green line algorithm
#     goal_size = len(network1) - 1

#     # init lists of sequence of moves
#     # list of (moving_edge,moving_endpoint,from_edge,to_edge)
#     seq_from_1 = []
#     seq_from_2 = []
#     # TODO keep track of lowest nodes?

#     # Do the green line algorithm
#     while (isom_size < goal_size):
#         # Find all lowest nodes above the isom in the networks:
#         lowest_tree_node_network1, lowest_retic_network1 = LowestReticAndTreeNodeAbove(network1, isom_1_2.keys(),
#                                                                                        allnodes=True)
#         lowest_tree_node_network2, lowest_retic_network2 = LowestReticAndTreeNodeAbove(network2, isom_2_1.keys(),
#                                                                                        allnodes=True)

#         # Construct a list of all lowest nodes in a tuple with the corresponding network (in random order)
#         # I.e. If u is a lowest node of network one, it will appear in the list as (u,1)
#         lowest_nodes_network1 = [(u, 1) for u in lowest_tree_node_network1 + lowest_retic_network1]
#         lowest_nodes_network2 = [(u, 2) for u in lowest_tree_node_network2 + lowest_retic_network2]
#         candidate_lowest_nodes = lowest_nodes_network1 + lowest_nodes_network2
#         random.shuffle(candidate_lowest_nodes)

#         # As some cases do not give an addition to the isom, we continue trying lowest nodes until we find one that does.
#         for lowest_node, network_number in candidate_lowest_nodes:
#             ######################################
#             # Case1: a lowest retic in network1
#             if network_number == 1 and network1.in_degree(lowest_node) == 2:
#                 # use notation as in the paper network1 = N', network2 = N, where ' is denoted p
#                 up = lowest_node
#                 if head_moves:
#                     moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_rSPR(
#                         network2, network1, up, isom_2_1, isom_1_2, randomNodes=True)
#                 else:
#                     moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case1_Tail(
#                         network2, network1, up, isom_2_1, isom_1_2, randomNodes=True)
#                     if added_node_network_1 == None:
#                         # The networks are non-isom networks with 2 leaves and 1 reticulation
#                         return False
#                 # This case always gives a node to add to the isom
#                 break

#             ######################################
#             # Case2: a lowest retic in network2
#             elif network_number == 2 and network2.in_degree(lowest_node) == 2:
#                 # use notation as in the paper network2 = N', network1 = N, where ' is denoted p
#                 up = lowest_node
#                 if head_moves:
#                     moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_rSPR(
#                         network1, network2, up, isom_1_2, isom_2_1, randomNodes=True)
#                 else:
#                     moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case1_Tail(
#                         network1, network2, up, isom_1_2, isom_2_1, randomNodes=True)
#                     if added_node_network_1 == None:
#                         # The networks are non-isom networks with 2 leaves and 1 reticulation
#                         return False
#                         # This case always gives a node to add to the isom
#                 break

#             ######################################
#             # Case3: a lowest tree node in network1
#             elif network_number == 2 and network2.out_degree(lowest_node) == 2:
#                 # use notation as in the paper network1 = N, network2 = N'
#                 up = lowest_node
#                 moves_network_1, moves_network_2, added_node_network_1, added_node_network_2 = GL_Case3(network1,
#                                                                                                         network2, up,
#                                                                                                         isom_1_2,
#                                                                                                         isom_2_1,
#                                                                                                         randomNodes=True)
#                 # If we can add a node to the isom, added_node_network_2 has a value
#                 if added_node_network_2:
#                     break

#             ######################################
#             # Case3': a lowest tree node in network2
#             else:
#                 # use notation as in the paper network1 = N, network2 = N'
#                 up = lowest_node
#                 moves_network_2, moves_network_1, added_node_network_2, added_node_network_1 = GL_Case3(network2,
#                                                                                                         network1, up,
#                                                                                                         isom_2_1,
#                                                                                                         isom_1_2,
#                                                                                                         randomNodes=True)
#                 # If we can add a node to the isom, added_node_network_2 has a value
#                 if added_node_network_2:
#                     break

#         # Now perform the moves and update the isomorphism
#         isom_1_2[added_node_network_1] = added_node_network_2
#         isom_2_1[added_node_network_2] = added_node_network_1
#         for move in moves_network_1:
#             seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
#             network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
#         for move in moves_network_2:
#             seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
#             network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
#         isom_size += 1
#     # TESTING FOR CORRECTNESS WHILE RUNNING
#     #        if not Isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
#     #            print("not unlabeled isom")
#     #            print(seq_from_1)
#     #            print(seq_from_2)
#     #            print(isom_1_2)
#     #            print(network1.subgraph(isom_1_2.keys()).edges())
#     #            print(network2.subgraph(isom_2_1.keys()).edges())

#     # Add the root to the isomorphism, if it was there
#     isom_1_2[root1] = root2
#     isom_2_1[root2] = root1

#     # invert seq_from_2, rename to node names of network1, and append to seq_from_1
#     return seq_from_1 + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)
