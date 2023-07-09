# import os
# import sys
# import math
# import re
# import ast
# import random
# import numpy as np
# import pandas as pd
# from pathlib import Path
# import multiprocessing
# from multiprocessing import Manager
# from ExtractNetworkProperties import *
# import itertools


# def Newick_To_Network(newick):
#     newick=newick[:-1]

#     #remove internal labels
#     newick = re.sub(r"I([\d]+)", "", newick)

#     #remove lengths
#     newick = re.sub(r"([\d]+)\.([\d]+)", "", newick)
#     newick = re.sub(r"E-[\d]+", "", newick)
#     newick = re.sub(r":", "", newick)

#     #make into list format
#     newick = newick.replace("(","[")
#     newick = newick.replace(")","]")
#     newick = re.sub(r"\]\#H([\d]+)", r",#R\1]", newick)
#     newick = re.sub(r"#([RH])([\d]+)", r"'#\1\2'", newick)


#     #add "" if necessary
#     newick = re.sub(r"([ABCD])", r"'\1'", newick)
#     newick = re.sub(r" ", "", newick)


#     nestedtree = ast.literal_eval(newick)
#     edges, leaves, label_set, current_node = NestedList_To_Tree(nestedtree,1)
#     edges.append([0,1])
#     ret_labels = dict()
#     leaf_labels = dict()
#     for l in leaves:
#         if len(l)>2 and (l[:2]=="#H" or l[:2]=="#R"):
#             ret_labels[l[2:]]=[]
#         else:
#             leaf_labels[l]=[]
#     for l in label_set:
#         if len(l[0])>2 and (l[0][:2]=="#H" or l[0][:2]=="#R"):
#             if l[0][1]=='H':
#                 ret_labels[l[0][2:]]+=[l[1]]
#             else:
#                 ret_labels[l[0][2:]]=[l[1]]+ret_labels[l[0][2:]]
#         else:
#             leaf_labels[l[0]]+=[l[1]]
#     network = nx.DiGraph()
#     network.add_edges_from(edges)
#     for retic in ret_labels:
#         r = ret_labels[retic]
#         receiving = r[0]
#         parent_receiving = 0
#         for p in network.predecessors(receiving):
#             parent_receiving = p
#         network.remove_node(receiving)
#         for v in r[1:]:
#             network.add_edge(v,parent_receiving)
#             network = nx.contracted_edge(network,(v,parent_receiving))
#             network.remove_edge(v,v)
#             parent_receiving = v
#     leaves = set()
#     for l in leaf_labels:
#          leaf_labels[l]=leaf_labels[l][0]
#          leaves.add(l)
#     return network, leaves, leaf_labels


# def NestedList_To_Tree(nestedList,next_node):
#     edges = []
#     leaves = set()
#     labels = []
#     top_node = next_node
#     current_node = next_node+1
#     for t in nestedList:
#         edges.append((top_node,current_node))
#         if type(t)==list:
#             extra_edges, extra_leaves, extra_labels, current_node = NestedList_To_Tree(t,current_node)
#         else:
#             extra_edges = []
#             extra_leaves = set([str(t)])
#             extra_labels = [[str(t), current_node]]
#             current_node+=1
#         edges = edges + extra_edges
#         leaves = leaves.union(extra_leaves)
#         labels = labels + extra_labels
#     return edges, leaves, labels, current_node


# network,_,leaf_labels = Newick_To_Network(newick_string)
# reverse_labels = {x:y for (y,x) in leaf_labels.items()}


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########                         I/O Functions                     #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# ########
# ######## Convert Newick to a networkx Digraph with labels (and branch lengths)
# ########
# # Write length newick: convert ":" to "," and then evaluate as list of lists using ast.literal_eval
# # Then, in each list, the node is followed by the length of the incoming arc.
# # This only works as long as each branch has length and all internal nodes are labeled.
# def Newick_To_Network(newick):
#     """
#     Converts a Newick string to a networkx DAG with leaf labels.

#     :param newick: a string in extended Newick format for phylogenetic networks.
#     :return: a phylogenetic network, i.e., a networkx digraph with leaf labels represented by the `label' node attribute.
#     """
#     # Ignore the ';'
#     newick = newick[:-1]
#     # If names are not in string format between ', add these.
#     if not "'" in newick and not '"' in newick:
#         newick = re.sub(r"\)#H([\d]+)", r",#R\1)", newick)
#         newick = re.sub(r"([,\(])([#a-zA-Z\d]+)", r"\1'\2", newick)
#         newick = re.sub(r"([#a-zA-Z\d])([,\(\)])", r"\1'\2", newick)
#     else:
#         newick = re.sub(r"\)#H([d]+)", r"'#R\1'\)", newick)
#     newick = newick.replace("(", "[")
#     newick = newick.replace(")", "]")
#     nestedtree = ast.literal_eval(newick)
#     edges, current_node = NestedList_To_Network(nestedtree, 0, 1)
#     network = nx.DiGraph()
#     network.add_edges_from(edges)
#     network = NetworkLeafToLabel(network)
#     return network


# # Returns a network in which the leaves have the original name as label, and all nodes have integer names.
# def NetworkLeafToLabel(network):
#     """
#     Renames the network nodes to integers, while storing the original node names in the `original' node attribute.

#     :param network: a phylogenetic network
#     :return: a phylogenetic network with original node names in the `original' node attribute.
#     """
#     for node in network.nodes():
#         if network.out_degree(node) == 0:
#             network.node[node]['label'] = node
#     return nx.convert_node_labels_to_integers(network, first_label=0, label_attribute='original')


# # Auxiliary function to convert list of lists to graph
# def NestedList_To_Network(nestedList, top_node, next_node):
#     """
#     Auxiliary function used to convert list of lists to graph.

#     :param nestedList: a nested list.
#     :param top_node: an integer, the node name of the top node of the network represented by the list.
#     :param next_node: an integer, the lowest integer not yet used as node name in the network.
#     :return: a set of edges of the network represented by the nested list, and an updated next_node.
#     """
#     edges = []
#     if type(nestedList) == list:
#         if type(nestedList[-1]) == str and len(nestedList[-1]) > 2 and nestedList[-1][:2] == '#R':
#             retic_node = '#H' + nestedList[-1][2:]
#             bottom_node = retic_node
#         else:
#             bottom_node = next_node
#             next_node += 1
#         edges.append((top_node, bottom_node))
#         for t in nestedList:
#             extra_edges, next_node = NestedList_To_Network(t, bottom_node, next_node)
#             edges = edges + extra_edges
#     else:
#         if not (len(nestedList) > 2 and nestedList[:2] == '#R'):
#             edges = [(top_node, nestedList)]
#     return edges, next_node


# # Sets the labels of the nodes of a network with a given label dictionary
# def Set_Labels(network, label_dict):
#     """
#     Sets the labels of the leaves of a network using a dictionary of labels.

#     :param network: a networkx digraph, a DAG.
#     :param label_dict: a dictionary, containing the labels (values) of the nodes of the network (keys).
#     :return: a phylogenetic network, obtained by labeling network with the labels.
#     """
#     for node, value in label_dict.items():
#         network.node[node]['label'] = value


################################################################################
################################################################################
################################################################################
########                                                           #############
########                     AAE CutTree CLASS                     #############
########                                                           #############
################################################################################
################################################################################
################################################################################


# #A class that represents a network as a tree where hybrid edges have been cut at the hybrid nodes.
# #Used as an intermediate to find the Newick string of a network.
# class CutTree:
#     def __init__(self, network = None, current_node = None, leaf_labels= dict()):
#          self.hybrid_nodes = dict()
#          self.no_of_hybrids = 0
#          self.root = None
#          self.nw = deepcopy(network)
#          self.current_node = current_node
#          self.leaf_labels = leaf_labels
#          if not self.current_node:
#              self.current_node = 2*len(self.nw)
#          if network:
#              self.Find_Root()
#              network_nodes = list(self.nw.nodes)
#              for node in network_nodes:
#                  if self.nw.in_degree(node)>1:
#                      self.no_of_hybrids+=1
#                      enumerated_parents = list(enumerate(self.nw.predecessors(node)))
#                      for i,parent in enumerated_parents:
#                          if i==0:
#                              self.hybrid_nodes[node]=self.no_of_hybrids
#                          else:
#                              self.nw.add_edges_from([(parent,self.current_node,self.nw[parent][node])])
#                              self.nw.remove_edge(parent,node)
#                              self.hybrid_nodes[self.current_node] = self.no_of_hybrids
#                              self.current_node+=1
# #             self.CheckLabelSet()

#     #Returns the root node of the tree
#     def Find_Root(self):
#         for node in self.nw.nodes:
#             if self.nw.in_degree(node)==0:
#                 self.root = node
#                 return node

#     #Returns a newick string for the tree
#     def Newick(self,probabilities = False):
#         return self.Newick_Recursive(self.root,probabilities = probabilities)+";"

#     #Returns the newick string for the subtree with given root
#     #does not append the; at the end, for the full newick string of the tree, use Newick()
#     # auxiliary function for finding the newick string for the tree
#     def Newick_Recursive(self,root,probabilities = False):
#         if self.nw.out_degree(root)==0:
#             if root in self.hybrid_nodes:
#                 return "#H"+str(self.hybrid_nodes[root])
#             elif root in self.leaf_labels:
#                 return self.leaf_labels[root]
#             return str(root)
#         Newick = ""
#         for v in self.nw.successors(root):
#             Newick+= self.Newick_Recursive(v,probabilities)+":"+str(self.nw[root][v]['length'])
#             if probabilities and v in self.hybrid_nodes:
#                 Newick+="::"+str(self.nw[root][v]['prob'])
#             Newick+= ","
#         Newick = "("+Newick[:-1]+")"
#         if root in self.hybrid_nodes:
#             Newick += "#H"+str(self.hybrid_nodes[root])
#         return Newick

#     '''
#     def CheckLabelSet(self):
#         for v in self.nw.nodes:
#             if self.nw.out_degree(v)==0:
#                 if v not in self.leaf_labels and v not in self.hybrid_nodes:
#                     print("non-labelled leaf!")
#                     return False
#         return True
#     '''
