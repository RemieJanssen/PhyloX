# import numpy as np
# import random
# import networkx as nx
# import sys


# ##PARAMETERS
# time_limit = 0.2
# speciation_rate = 20.0
# hybridization_rate = 1.0
# inheritence = True
# edges = True
# ##


# ###############################2. I/O############################

# i = 1
# while i < len(sys.argv):
#     arg= sys.argv[i]
#     if arg == "-t":
#         i+=1
#         time_limit = float(sys.argv[i])
#     if arg == "-sp":
#         i+=1
#         speciation_rate = float(sys.argv[i])
#     if arg == "-hyb":
#         i+=1
#         hybridization_rate = float(sys.argv[i])
#     if arg == "-no_inheritence":
#         inheritence = False
#     i += 1


# nw = nx.DiGraph()
# nw.add_node(0)
# leaves = set([0])
# current_node = 1

# extra_time = np.random.exponential(1/float(speciation_rate))
# current_time = extra_time
# current_speciation_rate    = float(speciation_rate)
# current_hybridization_rate = float(0)
# rate = current_speciation_rate + current_hybridization_rate

# #First create a MUL-tree
# hybrid_nodes=dict()
# no_of_hybrids = 0

# while current_time<time_limit:
#     if random.random() < current_speciation_rate / rate:
#         #speciate
#         splitting_leaf = random.choice(list(leaves))
#         nw.add_weighted_edges_from([(splitting_leaf,current_node,0),(splitting_leaf,current_node+1,0)], weight = 'length')
#         leaves.remove(splitting_leaf)
#         leaves.add(current_node)
#         leaves.add(current_node+1)
#         current_node+=2
#     else:
#         #Hybridize
# 	#i.e.: pick two leaf nodes, merge those, and add a new leaf below this hybrid node.
#         merging = random.sample(leaves,2)
#         l0 = merging[0]
#         l1 = merging[1]
#         pl0 = -1
#         for p in nw.predecessors(l0):
#             pl0=p
#         pl1 = -1
#         for p in nw.predecessors(l1):
#             pl1=p
#         #If pl0==pl1, the new hybridization results in parallel edges.
#         if pl0 != pl1:
#             no_of_hybrids+=1
#             nw.add_weighted_edges_from([(l0,current_node,0)],weight = 'length')
#             leaves.remove(l0)
#             leaves.remove(l1)
#             leaves.add(current_node)
#             prob = random.random()
#             nw[pl0][l0]['prob'] = prob
#             nw[pl1][l1]['prob'] = 1-prob
#             hybrid_nodes[l0]=no_of_hybrids
#             hybrid_nodes[l1]=no_of_hybrids
#             current_node+=1
#     #Now extend all pendant edges
#     for l in leaves:
#         pl = -1
#         for p in nw.predecessors(l):
#             pl = p
#         nw[pl][l]['length']+=extra_time
#     no_of_leaves = len(leaves)
#     current_speciation_rate    = float(speciation_rate*no_of_leaves)
#     current_hybridization_rate = float(hybridization_rate * (no_of_leaves * (no_of_leaves - 1))/2)
#     rate = current_speciation_rate + current_hybridization_rate
#     extra_time    = np.random.exponential(1/rate)
#     current_time += extra_time


# extra_time -= current_time-time_limit
# #nothing has happened yet, and there is only one node
# if len(nw) == 1:
#     nw.add_weighted_edges_from([(0,1,time_limit)],weight = 'length')
#     leaves = set([1])
# # each leaf has a parent node, and we can extend each parent edge to time_limit
# else:
#     for l in leaves:
#         pl = -1
#         for p in nw.predecessors(l):
#             pl = p
#         nw[pl][l]['length']+=extra_time


# ############### NOW CONVERT TO NEWICK ##############

# def Newick_From_MULTree(tree,root,hybrid_nodes):
#     if tree.out_degree(root)==0:
#         if root in hybrid_nodes:
#             return "#H"+str(hybrid_nodes[root])
#         return str(root)
#     Newick = ""
#     for v in tree.successors(root):
#         Newick+= Newick_From_MULTree(tree,v,hybrid_nodes)+":"+str(tree[root][v]['length'])
#         if inheritence:
#             if v in hybrid_nodes:
#                 Newick+="::"+str(tree[root][v]['prob'])
#         Newick+= ","
#     Newick = "("+Newick[:-1]+")"
#     if root in hybrid_nodes:
#         Newick += "#H"+str(hybrid_nodes[root])
#     return Newick


# print(Newick_From_MULTree(nw,0,hybrid_nodes)+";")
# inheritence = False
# print(Newick_From_MULTree(nw,0,hybrid_nodes)+";")
# if edges:
#     for e in nw.edges(data=True):
#         info = ""
#         if e[0] in hybrid_nodes:
#             info += "#H"+str(hybrid_nodes[e[0]])
#         else:
#             info += str(e[0])
#         info += ","
#         if e[1] in hybrid_nodes:
#             info += "#H"+str(hybrid_nodes[e[1]])
#         else:
#             info += str(e[1])
#         info += ","
#         info += str(e[2]['length']) + ",1.0"
#         if 'prob' in e[2]:
#             info = info[:-3] + str(e[2]['prob'])
#         print(info)
