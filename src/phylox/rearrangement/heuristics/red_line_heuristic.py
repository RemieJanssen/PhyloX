# # Code used for cases 1 and 2 in the head move red line algorithm
# # returns moves1,moves2,addedNode1,addedNode2
# def RL_Case1(N1, N2, x_1, isom_N1_N2, isom_N2_N1, randomNodes=False):
#     """
#     An implementation of Algorithm 7. Finds a sequence of head moves that makes it possible to add the highest tree node x_1 to the up-closed isomrophism.

#     :param N1: a phylogenetic network.
#     :param N2: a phylogenetic network.
#     :param x_1: a tree node in N1.
#     :param isom_N1_N2: a dictionary, containing a partial (up-closed) isomorphism map from N1 to N2. The inverse of isom_N2_N1.
#     :param isom_N2_N1: a dictionary, containing a partial (up-closed) isomorphism map from N2 to N1. The inverse of isom_N1_N2.
#     :param randomNodes: a boolean value, determining whether the random version of this algorithm is used.
#     :return: a list of head moves in N1, a list of head moves in N2, a node of N1, a node of N2. After applying the moves to the networks, the nodes can be added to the up-closed isomorphism.
#     """
#     p_1 = Parent(N1, x_1)
#     p_2 = isom_N1_N2[p_1]
#     x_2 = Child(N2, p_2, exclude=isom_N2_N1.keys(), randomNodes=randomNodes)
#     if N2.out_degree(x_2) == 2:
#         #        print("Case tree")
#         return [], [], x_1, x_2
#     elif N2.in_degree(x_2) == 2:
#         #        print("Case retic")
#         c_2 = FindTreeNode(N2, excludedSet=isom_N2_N1.keys(), randomNodes=randomNodes)
#         t_2 = Parent(N2, c_2)
#         b_2 = Child(N2, c_2, exclude=[x_2], randomNodes=randomNodes)
#         #        print(x_1,p_1,p_2,x_2,c_2,t_2,b_2)
#         # TODO: Not in the latex algo, just a minor improvement:
#         # If the other child $c_2$ of $p_2$ is a retic, then we can add it.
#         if p_2 == t_2:
#             return [], [], x_1, c_2
#         if check_movable(N2, (p_2, x_2), x_2):
#             q_2 = Parent(N2, x_2, exclude=[p_2])
#             if x_2 == t_2:
#                 #                print("(p,x) movable and x=t")
#                 return [], [((q_2, x_2), x_2, (c_2, b_2))], x_1, c_2
#             else:
#                 #                print("(p,x) movable and x!=t")
#                 return [], [((p_2, x_2), x_2, (t_2, c_2)), ((t_2, x_2), x_2, (c_2, b_2))], x_1, c_2
#         else:
#             d_2 = Child(N2, x_2)
#             z_2 = Parent(N2, x_2, exclude=[p_2])
#             # Find a leaf with parent not equal to d_2
#             l_2 = FindLeaf(N2, excludedParents=[d_2], randomNodes=randomNodes)
#             w_2 = Parent(N2, l_2)
#             #            print("other")
#             #            print(d_2,z_2,w_2,l_2)
#             if l_2 == b_2:  # after the first move, b_2 is not a child of c_2 anymore, so we fix this by taking the other child of c_2 as b_2
#                 b_2 = Child(N2, c_2, exclude=[l_2, x_2], randomNodes=randomNodes)
#             if d_2 == t_2:  # i.e., if, after the first move, x_2 is the parent of c_2
#                 return [], [((z_2, d_2), d_2, (w_2, l_2)), ((z_2, x_2), x_2, (c_2, b_2))], x_1, c_2
#             else:
#                 return [], [((z_2, d_2), d_2, (w_2, l_2)), ((p_2, x_2), x_2, (t_2, c_2)),
#                             ((t_2, x_2), x_2, (c_2, b_2))], x_1, c_2
#     else:
#         #        print("Case leaf")
#         c_2 = FindTreeNode(N2, excludedSet=isom_N2_N1.keys(), randomNodes=randomNodes)
#         t_2 = Parent(N2, c_2)

#         if p_2 == t_2:
#             #            print("no move")
#             return [], [], x_1, c_2

#             # Find a reticulation arc (s_2,r_2) that can be moved to (p_2,x_2)
#         # No randomness required, because this arc will end up at its original position again.
#         s_2 = None
#         r_2 = None
#         for node in N2.nodes():
#             if N2.in_degree(node) == 2:
#                 for parent in N2.predecessors(node):
#                     if parent != p_2 and check_movable(N2, (parent, node), node):
#                         s_2 = parent
#                         r_2 = node
#             if s_2 != None:
#                 break
#         q_2 = Parent(N2, r_2, exclude=[s_2])
#         w_2 = Child(N2, r_2)
#         #        print(p_2,x_2,t_2,c_2,s_2,q_2,r_2,w_2)
#         if r_2 == p_2:
#             #            print("r=p")
#             if s_2 == t_2:
#                 return [], [((q_2, p_2), p_2, (t_2, c_2))], x_1, c_2
#             elif q_2 == t_2:
#                 return [], [((s_2, p_2), p_2, (t_2, c_2))], x_1, c_2
#             else:
#                 return [], [((s_2, p_2), p_2, (t_2, c_2)), ((t_2, p_2), p_2, (q_2, x_2)),
#                             ((q_2, p_2), p_2, (s_2, c_2))], x_1, c_2
#         if r_2 == t_2:
#             #            print("r=t")
#             if p_2 == q_2:
#                 return [], [((s_2, r_2), r_2, (p_2, x_2))], x_1, c_2
#             else:
#                 return [], [((s_2, r_2), r_2, (p_2, x_2)), ((p_2, r_2), r_2, (q_2, c_2)),
#                             ((q_2, r_2), r_2, (s_2, x_2))], x_1, c_2
#         if s_2 != t_2:
#             #            print("s!=t")
#             return [], [((s_2, r_2), r_2, (p_2, x_2)), ((p_2, r_2), r_2, (t_2, c_2)), ((t_2, r_2), r_2, (s_2, x_2)),
#                         ((s_2, r_2), r_2, (q_2, w_2))], x_1, c_2
#         else:
#             #            print("other")
#             return [], [((s_2, r_2), r_2, (p_2, x_2)), ((p_2, r_2), r_2, (s_2, c_2)),
#                         ((s_2, r_2), r_2, (q_2, w_2))], x_1, c_2

#         # Code used for case 3 in the head move red line algorithm


# # returns moves1,moves2,addedNode1,addedNode2
# def RL_Case3(N1, N2, x_1, isom_N1_N2, isom_N2_N1, randomNodes=False):
#     """
#     An implementation of Algorithm 8. Finds a sequence of head moves that makes it possible to add the highest reticulation x_1 to the up-closed isomrophism.

#     :param N1: a phylogenetic network.
#     :param N2: a phylogenetic network.
#     :param x_1: a highest reticulation node in N1.
#     :param isom_N1_N2: a dictionary, containing a partial (up-closed) isomorphism map from N1 to N2. The inverse of isom_N2_N1.
#     :param isom_N2_N1: a dictionary, containing a partial (up-closed) isomorphism map from N2 to N1. The inverse of isom_N1_N2.
#     :param randomNodes: a boolean value, determining whether the random version of this algorithm is used.
#     :return: a list of head moves in N1, a list of head moves in N2, a node of N1, a node of N2. After applying the moves to the networks, the nodes can be added to the up-closed isomorphism.
#     """
#     p_1 = Parent(N1, x_1)
#     q_1 = Parent(N1, x_1, exclude=[p_1])
#     p_2 = isom_N1_N2[p_1]
#     cp_2 = Child(N2, p_2, exclude=isom_N2_N1.keys(), randomNodes=randomNodes)
#     q_2 = isom_N1_N2[q_1]
#     cq_2 = Child(N2, q_2, exclude=isom_N2_N1.keys(), randomNodes=randomNodes)

#     # at least one tree node, so give up
#     if N2.out_degree(cp_2) == 2 or N2.out_degree(cq_2) == 2:
#         # The proof does not provide a sequence when at least of the nodes cp_2 or cq_2 is a tree node
#         # TODO: If one of (p_2,cp_2) or (q_2,cq_2) is movable, we can still do something quite similar to what follows in the last case
#         return [], [], None, None
#     # Case 3ai
#     if cp_2 == cq_2:
#         #        print("no move")
#         return [], [], x_1, cp_2
#     # Case 3av
#     elif N2.out_degree(cp_2) == 0 and N2.out_degree(cq_2) == 0:
#         #        print("both leaves")
#         # Find a head-movable arc (s_2,r_2)
#         # As each reticulation node has a movable arc, we pick a random reticulation, and then a random movable arc incident to this node
#         r_2 = FindRetic(N2, excludedSet=isom_N2_N1.keys(), randomNodes=randomNodes)
#         s_2 = None
#         for parent in N2.predecessors(r_2):
#             if check_movable(N2, (parent, r_2), r_2):
#                 if not randomNodes:
#                     s_2 = parent
#                     break
#                 elif s_2 == None or random.getrandbits(1):
#                     s_2 = parent
#         if s_2 == p_2:
#             #            print("s=p")
#             return [], [((s_2, r_2), r_2, (q_2, cq_2)), ((q_2, r_2), r_2, (p_2, cp_2))], x_1, r_2
#         else:
#             #            print("s!=p")
#             return [], [((s_2, r_2), r_2, (p_2, cp_2)), ((p_2, r_2), r_2, (q_2, cq_2))], x_1, r_2
#     # Cases 3a(ii, iii, iv)
#     else:
#         #        print("at least one non-leaf")
#         if nx.has_path(N2, cq_2, cp_2) or N2.out_degree(cp_2) == 0:
#             # Swap p and q
#             q_2, p_2 = p_2, q_2
#             cp_2, cq_2 = cq_2, cp_2
#         #        print(p_2 , q_2, cp_2, cq_2)
#         if check_movable(N2, (p_2, cp_2), cp_2):
#             #            print("movable cp_2")
#             return [], [((p_2, cp_2), cp_2, (q_2, cq_2))], x_1, cp_2
#         else:
#             #            print("immovable cp_2")
#             z = Child(N2, cp_2)
#             t = Parent(N2, cp_2, exclude=[p_2])
#             if t == q_2:
#                 return [], [], x_1, cp_2
#             return [], [((cp_2, z), z, (q_2, cq_2)), ((t, cp_2), cp_2, (z, cq_2))], x_1, z


# # Permutes the leaves of a network so that it becomes leaf isomorphic, provided the networks were already (non-labeled) ismorphic
# def Permute_Leaves_Head(network1, network2, isom_1_2, isom_2_1, label_dict_1, label_dict_2):
#     """
#     An implementation of Algorithm 9. Determines a sequence of head moves that makes two isomorphic networks labeled isomorphic.

#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :param isom_1_2: a dictionary, containing a partial (up-closed) isomorphism map from network1 to network2. The inverse of isom_2_1.
#     :param isom_2_1: a dictionary, containing a partial (up-closed) isomorphism map from network2 to network1. The inverse of isom_1_2.
#     :param label_dict_1: a dictionary, containing the correspondence of nodes of network1 with the leaf labels: keys are labels and values are nodes.
#     :param label_dict_2: a dictionary, containing the correspondence of nodes of network2 with the leaf labels: keys are labels and values are nodes.
#     :return: a list of head moves that turns network1 into network2.
#     """
#     #    for i,k in isom_1_2.items():
#     #        print(i,k)
#     sequence = []
#     # labeldict[label]=leaf
#     Y = list(label_dict_1.values())
#     cycles = []
#     while len(Y) > 0:
#         y1_1 = Y.pop()
#         y_2 = isom_1_2[y1_1]
#         cycle = [y1_1]
#         while network2.node[y_2].get('label') != network1.node[cycle[0]].get('label'):
#             y_new_1 = label_dict_1[network2.node[y_2].get('label')]
#             #            if len(set(network1.predecessors(cycle[-1]))&set(network1.predecessors(y_new_1)))==0:#cycle[-1] and y_new_1 have NO common parent
#             #                cycle+=[y_new_1]
#             cycle += [y_new_1]

#             y_2 = isom_1_2[y_new_1]
#             Y.remove(y_new_1)
#         if len(cycle) > 1:
#             cycles += [list(reversed(cycle))]
#     #    print("cycles",cycles)

#     t = None
#     r = None
#     for node in network1:
#         if network1.in_degree(node) == 2:
#             for parent in network1.predecessors(node):
#                 if check_movable(network1, (parent, node), node):
#                     t = parent
#                     r = node
#         if r:
#             break
#     c_last = Child(network1, r)  # In the proof: c', the child of r in N
#     s_last = Parent(network1, r, exclude=[t])

#     for cycle in cycles:
#         c = cycle

#         # shift the cycle by one of t is the parent of the last leaf in the cycle
#         if t in network1.predecessors(cycle[-1]):
#             c = [cycle[-1]] + cycle[:-1]

#         # Skip the first move if the head r of the moving arc is already above the last leaf in the cycle
#         p_last = Parent(network1, c[-1])
#         if r != p_last:
#             move = ((t, r), r, (p_last, c[-1]))
#             sequence.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
#             network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
#         moved_arc = (t, r)

#         c_last_before = c_last

#         for i in reversed(range(len(c))):
#             p_i = Parent(network1, r, exclude=[moved_arc[0]])
#             p_im1 = Parent(network1, c[i - 1])
#             # if p_i==p_im1, then two consecutive leaves are in a cherry, and we can skip one of them.
#             if p_i == p_im1:
#                 continue
#             #                print("do nothing, swapping a cherry")
#             else:
#                 move = ((p_i, r), r, (p_im1, c[i - 1]))
#                 sequence.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
#                 network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
#                 moved_arc = (p_i, r)
#             # If c in the `original position' (s,c) of (t,r) is permuted by the cycle, then we need to change this original position accordingly
#             if c[i] == c_last_before:
#                 c_last = c[i - 1]

#     # Put the moving arc back to its original position
#     # The proof does this for every cycle, but can also be done once at the end
#     if (s_last, c_last) in network1.edges():
#         # (t,r) might already be back at this position
#         move = ((t, r), r, (s_last, c_last))
#         sequence.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
#         network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)

#     #    print("isomorphic", Isomorphic(network1,network2))
#     #    for s in sequence:
#     #        print(s[0])
#     return sequence


# # A heuristic for finding a head move sequence between two networks.
# def Red_Line(network1, network2):
#     """
#     An implementation of Algorithm 10. Finds a sequence of head moves from network1 to network2 by building an up-closed isomorphism.

#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :return: A list of head moves from network1 to network2.
#     """
#     # Find the root and labels of the networks
#     root1 = Root(network1)
#     root2 = Root(network2)
#     label_dict_1 = Labels(network1)
#     label_dict_2 = Labels(network2)

#     # initialize isomorphism
#     isom_1_2 = dict()
#     isom_1_2[root1] = root2
#     isom_2_1 = dict()
#     isom_2_1[root2] = root1
#     isom_size = 1

#     # Check if the roots are of the same type
#     if network1.out_degree(root1) != network2.out_degree(root2):
#         return False

#     # Keep track of the size of the isomorphism and the size it is at the end of the red line algorithm
#     goal_size = len(network1) - len(label_dict_1)

#     # init lists of sequence of moves
#     # list of (moving_edge,moving_endpoint,from_edge,to_edge)
#     seq_from_1 = []
#     seq_from_2 = []
#     # TODO keep track of highest nodes?

#     # Do the red line algorithm
#     while (isom_size < goal_size):
#         highest_tree_node_network1, highest_retic_network1 = HighestNodesBelow(network1, isom_1_2.keys())
#         highest_tree_node_network2, highest_retic_network2 = HighestNodesBelow(network2, isom_2_1.keys())

#         # Case1
#         if highest_tree_node_network1 != None:
#             moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case1(network1, network2,
#                                                                                                   highest_tree_node_network1,
#                                                                                                   isom_1_2, isom_2_1)
#         # Case2
#         elif highest_tree_node_network2 != None:
#             moves_network_2, moves_network_1, added_node_network2, added_node_network1 = RL_Case1(network2, network1,
#                                                                                                   highest_tree_node_network2,
#                                                                                                   isom_2_1, isom_1_2)
#         # Case3
#         else:
#             moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case3(network1, network2,
#                                                                                                   highest_retic_network1,
#                                                                                                   isom_1_2, isom_2_1)

#         # Now perform the moves and update the isomorphism
#         isom_1_2[added_node_network1] = added_node_network2
#         isom_2_1[added_node_network2] = added_node_network1
#         for move in moves_network_1:
#             seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
#             network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
#         for move in moves_network_2:
#             seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
#             network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
#         isom_size += 1
#     # TESTING FOR CORRECTNESS WHILE RUNNING
#     #        if not nx.is_isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
#     #            print("not unlabeled isom")
#     #            print(seq_from_1)
#     #            print(seq_from_2)
#     #            print(network1.subgraph(isom_1_2.keys()).edges())
#     #            print(network2.subgraph(isom_2_1.keys()).edges())

#     # TODO Debugging, remove after for speed
#     #    if not nx.is_isomorphic(network1,network2):
#     #        print("not unlabeled isom")
#     #        print(network1.edges())
#     #        print(network2.edges())
#     #    else:
#     #        print("unlabeled isomorphic :)")

#     # Add the leaves to the isomorphism
#     for node_1 in network1.nodes():
#         if network1.out_degree(node_1) == 0:
#             parent_1 = Parent(network1, node_1)
#             parent_2 = isom_1_2[parent_1]
#             node_2 = Child(network2, parent_2, exclude=isom_2_1.keys())
#             isom_1_2[node_1] = node_2
#             isom_2_1[node_2] = node_1

#     # Permute the leaves
#     seq_permute = Permute_Leaves_Head(network1, network2, isom_1_2, isom_2_1, label_dict_1, label_dict_2)

#     # invert seq_from_2, rename to node names of network1, and append to seq_from_1
#     return seq_from_1 + seq_permute + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)


# def Red_Line_Random(network1, network2, repeats=1):
#     """
#     Finds a sequence of head moves from network1 to network2 by randomly building an up-closed isomorphism a number of times, and only keeping the shortest sequence.

#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :param repeats: an integer, determining how many repeats of Red_Line_Random_Single are performed.
#     :return: A list of head moves from network1 to network2.
#     """
#     best_seq = None
#     for i in range(repeats):
#         candidate_seq = Red_Line_Random_Single(network1, network2)
#         if best_seq == None or len(best_seq) > len(candidate_seq):
#             best_seq = candidate_seq
#     return best_seq


# def Red_Line_Random_Single(network1, network2):
#     """
#     An implementation of Algorithm 11. Finds a sequence of head moves from network1 to network2 by building randomly an up-closed isomorphism.

#     :param network1: a phylogenetic network.
#     :param network2: a phylogenetic network.
#     :return: A list of head moves from network1 to network2.
#     """
#     # Find the root and labels of the networks
#     root1 = Root(network1)
#     root2 = Root(network2)
#     label_dict_1 = Labels(network1)
#     label_dict_2 = Labels(network2)

#     # initialize isomorphism
#     isom_1_2 = dict()
#     isom_1_2[root1] = root2
#     isom_2_1 = dict()
#     isom_2_1[root2] = root1
#     isom_size = 1

#     # Check if the roots are of the same type
#     if network1.out_degree(root1) != network2.out_degree(root2):
#         return False

#     # Keep track of the size of the isomorphism and the size it is at the end of the red line algorithm
#     goal_size = len(network1) - len(label_dict_1)

#     # init lists of sequence of moves
#     # list of (moving_edge,moving_endpoint,from_edge,to_edge)
#     seq_from_1 = []
#     seq_from_2 = []
#     # TODO keep track of highest nodes?

#     # Do the red line algorithm
#     while (isom_size < goal_size):
#         highest_tree_node_network1, highest_retic_network1 = HighestNodesBelow(network1, isom_1_2.keys(), allnodes=True)
#         highest_tree_node_network2, highest_retic_network2 = HighestNodesBelow(network2, isom_2_1.keys(), allnodes=True)

#         # Construct a list of all highest nodes in a tuple with the corresponding network (in random order)
#         # I.e. If u is a highest node of network one, it will appear in the list as (u,1)
#         highest_nodes_network1 = [(u, 1) for u in highest_tree_node_network1 + highest_retic_network1]
#         highest_nodes_network2 = [(u, 2) for u in highest_tree_node_network2 + highest_retic_network2]
#         candidate_highest_nodes = highest_nodes_network1 + highest_nodes_network2
#         random.shuffle(candidate_highest_nodes)

#         # As some cases do not give an addition to the isom, we continue trying lowest nodes until we find one that does.
#         for highest_node, network_number in candidate_highest_nodes:
#             u = highest_node
#             # Case1
#             if network_number == 1 and network1.in_degree(u) == 1:
#                 #                print("Case 1")
#                 moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case1(network1,
#                                                                                                       network2, u,
#                                                                                                       isom_1_2,
#                                                                                                       isom_2_1,
#                                                                                                       randomNodes=True)
#                 break
#             # Case2
#             elif network_number == 2 and network2.in_degree(u) == 1:
#                 #                print("Case 2")
#                 moves_network_2, moves_network_1, added_node_network2, added_node_network1 = RL_Case1(network2,
#                                                                                                       network1, u,
#                                                                                                       isom_2_1,
#                                                                                                       isom_1_2,
#                                                                                                       randomNodes=True)
#                 break
#             # Case3
#             elif network_number == 1 and network1.in_degree(u) == 2:
#                 #                print("Case 3")
#                 moves_network_1, moves_network_2, added_node_network1, added_node_network2 = RL_Case3(network1,
#                                                                                                       network2, u,
#                                                                                                       isom_1_2,
#                                                                                                       isom_2_1,
#                                                                                                       randomNodes=True)
#                 if added_node_network2:
#                     break
#             # Case3'
#             elif network_number == 2 and network2.in_degree(u) == 2:
#                 #                print("Case 3'")
#                 moves_network_2, moves_network_1, added_node_network2, added_node_network1 = RL_Case3(network2,
#                                                                                                       network1, u,
#                                                                                                       isom_2_1,
#                                                                                                       isom_1_2,
#                                                                                                       randomNodes=True)
#                 if added_node_network2:
#                     break

#         # Now perform the moves and update the isomorphism
#         isom_1_2[added_node_network1] = added_node_network2
#         isom_2_1[added_node_network2] = added_node_network1
#         for move in moves_network_1:
#             seq_from_1.append((move[0], move[1], From_Edge(network1, move[0], move[1]), move[2]))
#             network1 = DoMove(network1, move[0], move[1], move[2], check_validity=False)
#         for move in moves_network_2:
#             seq_from_2.append((move[0], move[1], From_Edge(network2, move[0], move[1]), move[2]))
#             network2 = DoMove(network2, move[0], move[1], move[2], check_validity=False)
#         isom_size += 1
#     # TESTING FOR CORRECTNESS WHILE RUNNING
#     #        if not nx.is_isomorphic(network1.subgraph(isom_1_2.keys()),network2.subgraph(isom_2_1.keys())):
#     #            print("not unlabeled isom")
#     #            print(seq_from_1)
#     #            print(seq_from_2)
#     #            print(network1.subgraph(isom_1_2.keys()).edges())
#     #            print(network2.subgraph(isom_2_1.keys()).edges())
#     #            print(network1.edges())
#     #            print(network2.edges())

#     # TODO Debugging, remove after for speed
#     #    if not nx.is_isomorphic(network1,network2):
#     #        print("not unlabeled isom")
#     #        print(network1.edges())
#     #        print(network2.edges())
#     #    else:
#     #        print("unlabeled isomorphic :)")

#     # Add the leaves to the isomorphism
#     for node_1 in network1.nodes():
#         if network1.out_degree(node_1) == 0:
#             parent_1 = Parent(network1, node_1)
#             parent_2 = isom_1_2[parent_1]
#             node_2 = Child(network2, parent_2, exclude=isom_2_1.keys())
#             isom_1_2[node_1] = node_2
#             isom_2_1[node_2] = node_1

#     # Permute the leaves
#     seq_permute = Permute_Leaves_Head(network1, network2, isom_1_2, isom_2_1, label_dict_1, label_dict_2)

#     # invert seq_from_2, rename to node names of network1, and append to seq_from_1
#     return seq_from_1 + seq_permute + ReplaceNodeNamesInMoveSequence(InvertMoveSequence(seq_from_2), isom_2_1)
