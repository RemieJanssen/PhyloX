# import networkx as nx
# import random
# import matplotlib.pyplot as plt
# import os
# import ast
# import re
# import time

# random.seed(a=1234321)

# ### For the tree child network containment paper by Murakami and Janssen

# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########                  PHYLOGENETIC NETWORK CLASS               #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# # A class for phylogenetic networks
# class PhN:
#     def __init__(self, seq=None, newick=None):
#         # the actual graph
#         self.nw = nx.DiGraph()
#         # the set of leaf labels of the network
#         self.leaves = set()
#         # a dictionary giving the node for a given leaf label
#         self.labels = dict()
#         # the number of nodes in the graph
#         self.no_nodes = 0
#         self.leaf_nodes = dict()
#         self.TCS = seq
#         self.CPS = seq
#         self.newick = newick
#         self.reducible_pairs = set()
#         self.reticulated_cherries = set()
#         self.cherries = set()
#         if seq:
#             # Creates a phylogenetic network from a cherry picking sequence:
#             for pair in reversed(seq):
#                 self.add_pair(*pair)

#     def Compute_Leaf_Nodes(self):
#         self.leaf_nodes = dict()
#         for v in self.labels:
#             self.leaf_nodes[self.labels[v]] = v

#     # A method for adding a pair, using the construction from a sequence
#     def add_pair(self, x, y):
#         if len(self.leaves) == 0:
#             self.nw.add_edges_from([(0, 1), (1, 2), (1, 3)])
#             self.leaves = set([x, y])
#             self.labels[x] = 2
#             self.labels[y] = 3
#             self.leaf_nodes[2] = x
#             self.leaf_nodes[3] = y
#             self.no_nodes = 4
#             return True
#         if y not in self.leaves:
#             return False
#         node_y = self.labels[y]
#         if x not in self.leaves:
#             self.nw.add_edges_from([(node_y, self.no_nodes), (node_y, self.no_nodes + 1)])
#             self.leaves.add(x)
#             self.leaf_nodes.pop(self.labels[y], False)
#             self.labels[y] = self.no_nodes
#             self.labels[x] = self.no_nodes + 1
#             self.leaf_nodes[self.no_nodes] = y
#             self.leaf_nodes[self.no_nodes + 1] = x
#             self.no_nodes += 2
#         else:
#             node_x = self.labels[x]
#             for parent in self.nw.predecessors(node_x):
#                 px = parent
#             if self.nw.in_degree(px) > 1:
#                 self.nw.add_edges_from([(node_y, px), (node_y, self.no_nodes)])
#                 self.leaf_nodes.pop(self.labels[y], False)
#                 self.labels[y] = self.no_nodes
#                 self.leaf_nodes[self.no_nodes] = y
#                 self.no_nodes += 1
#             else:
#                 self.nw.add_edges_from([(node_y, node_x), (node_y, self.no_nodes), (node_x, self.no_nodes + 1)])
#                 self.leaf_nodes.pop(self.labels[x], False)
#                 self.leaf_nodes.pop(self.labels[y], False)
#                 self.labels[y] = self.no_nodes
#                 self.labels[x] = self.no_nodes + 1
#                 self.leaf_nodes[self.no_nodes] = y
#                 self.leaf_nodes[self.no_nodes + 1] = x
#                 self.no_nodes += 2
#         return True


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########                       RANDOM NETWORKS                     #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# # A function that returns a tree-child sequence with a given number of leaves and reticulations
# def random_TC_sequence(leaves, retics):
#     current_leaves = set([1, 2])
#     seq = [(2, 1)]
#     not_forbidden = set([2])
#     leaves_left = leaves - 2
#     retics_left = retics

#     # Continue untill we added enough leaves and reticulations
#     while leaves_left > 0 or retics_left > 0:
#         # Decide if we add a leaf, or a reticulation
#         type_added = 'L'
#         if len(not_forbidden) > 0 and leaves_left > 0 and retics_left > 0:
#             if random.randint(0,
#                               leaves_left + retics_left - 1) < retics_left:  # probability of retic depends on number of retics left to add
#                 #           if random.randint(0 , 1)<1:                                        #probability of retics and leaves are the same if both are an option
#                 type_added = 'R'
#         elif len(not_forbidden) > 0 and retics_left > 0:
#             type_added = 'R'
#         elif leaves_left > 0:
#             type_added = 'L'
#         else:
#             return (False)

#         # Actually add the pair
#         if type_added == 'R':
#             first_element = random.choice(list(not_forbidden))
#             retics_left -= 1
#         if type_added == 'L':
#             first_element = len(current_leaves) + 1
#             leaves_left -= 1
#             current_leaves.add(first_element)
#             not_forbidden.add(first_element)

#         second_element = random.choice(list(current_leaves - set([first_element])))
#         not_forbidden.discard(second_element)
#         seq.append((first_element, second_element))

#     # reverse the sequence, as it was built in reverse order
#     seq = [pair for pair in reversed(seq)]
#     return (seq)


# # A function that returns a tree-child subsequence with a given number of reticulations
# def random_TC_subsequence(seq, r):
#     # First `uniformly at random' choose one pair per leaf, with that leaf as first element
#     leaves = dict()
#     indices = set()
#     for i, pair in enumerate(seq):
#         x = pair[0]
#         if x not in leaves:
#             indices.add(i)
#             leaves[x] = (1, i)
#         else:
#             if random.randint(0, leaves[x][0]) < 1:
#                 indices.remove(leaves[x][1])
#                 indices.add(i)
#                 leaves[x] = (leaves[x][0] + 1, i)
#             else:
#                 leaves[x] = (leaves[x][0] + 1, leaves[x][1])
#     # Add r reticulations with a max of the whole sequence
#     unused = set(range(len(seq))) - indices
#     for j in range(r):
#         new = random.choice(list(unused))
#         unused = unused - set([new])
#         indices.add(new)
#     newSeq = []
#     for i, pair in enumerate(seq):
#         if i in indices:
#             newSeq.append(pair)
#     return newSeq


# def make_a_lot_of_random_files(folder_name, edges=False):
#     try:
#         os.mkdir("./" + folder_name + "/")
#     except:
#         pass
#     f = open("./" + folder_name + "/hits_and_misses.txt", "w+")
#     f.write("index; subnetwork? \n")
#     f.close()
#     start = time.time()

#     maxLeavesRetics = 1000
#     stepLeavesRetics = 25
#     reps = 4

#     for leaves in range(stepLeavesRetics, maxLeavesRetics+1,stepLeavesRetics):
#         for reticulations in range(stepLeavesRetics, maxLeavesRetics+1,stepLeavesRetics):
#             for reticulationsSubnetwork in range(stepLeavesRetics, reticulations+1,stepLeavesRetics):
#                 for rep in range(reps):

#                     index1 = "0000000" + str(leaves)
#                     index1 = index1[-4:]
#                     index2 = "0000000" + str(reticulations)
#                     index2 = index2[-4:]
#                     index3 = "0000000" + str(reticulationsSubnetwork)
#                     index3 = index3[-4:]
#                     index = "n="+index1+"_k=" + index2+ "_kSub="+index3 + "_rep=" + str(rep)

#                     print(index)
#                     if os.path.isfile("./" + folder_name + "/" + index + ".txt"):
#                         print("already exists")
#                     else:
#                         sequence = random_TC_sequence(leaves, reticulations)
#                         snw_or_not = 'yes'
#                         if rep < reps/2:
#                             subsequence = random_TC_sequence(leaves, reticulations)
#                             snw_or_not = 'no'
#                         else:
#                             subsequence = random_TC_subsequence(sequence, reticulationsSubnetwork)
#                         f = open("./" + folder_name + "/" + index + ".txt", "w+")
#                         if edges:
#                             # build network from sequences
#                             network = PhN(seq=sequence)
#                             subnetwork = PhN(seq=subsequence)
#                             # now change the labels
#                             for node in network.leaf_nodes:
#                                 network.leaf_nodes[node] = "L" + str(network.leaf_nodes[node])
#                             for node in subnetwork.leaf_nodes:
#                                 subnetwork.leaf_nodes[node] = "L" + str(subnetwork.leaf_nodes[node])
#                             network.nw = nx.relabel_nodes(network.nw, network.leaf_nodes)
#                             subnetwork.nw = nx.relabel_nodes(subnetwork.nw, subnetwork.leaf_nodes)
#                             # Write edges to file
#                             # f.write(str(network.nw.edges)+"\r\n")
#                             f.write(str(network.nw.edges()) + "\r\n")
#                             f.write(str(subnetwork.nw.edges()))
#                         else:
#                             # Write newick to file
#                             f.write(SeqToNewick(sequence) + "\r\n")
#                             f.write(SeqToNewick(subsequence))
#                         f.close()
#                         # write answer to seperate file
#                         f = open("./" + folder_name + "/hits_and_misses.txt", "a+")
#                         f.write(index + " ; " + snw_or_not + "\n")
#                         f.close()
#     end = time.time()
#     print("time elapsed:", end - start, "seconds")


# #Name of the folder to put all input files in
# foldername = "tests"
# #make all the random input files
# make_a_lot_of_random_files(foldername, edges=True)
