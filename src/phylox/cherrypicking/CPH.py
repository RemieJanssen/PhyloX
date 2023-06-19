# import networkx as nx
# import random
# import matplotlib.pyplot as plt
# import os
# import ast
# import re
# import sys
# import time
# from copy import deepcopy
# memodict={}
# import numpy as np


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########                             Index                         #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################

# '''
# AAA I/O Functions
# 	Functions for reading newick sequences
# AAB ANALYSIS OF SEQUENCE FOR TREES
# 	Two functions:
#           - Checking whether a CPS reduces a set of trees
#           - Sequence_Add_Roots (i.e. CompletePartialSequence from the paper) to extend a partial CPS to a CPS
# AAC INPUT SET CLASS with CPS methods
# 	Class containing a set of inputs
#         Methods for running the CP heuristic and improving the sequence
# AAD PHYLOGENETIC TREE CLASS
# 	Class for a phylogenetic tree
# 	Contains methods to cherry-pick a tree and to find reducible pairs
# 	Gives the height of a pair in the tree
# AAE CutTree CLASS
# 	Class meant for converting a network into a Newick string
#         It `cuts' the reticulation arcs to produce a tree that can be converted to Newick
# AAF PHYLOGENETIC NETWORK CLASS
#         Class for a phylogenetic network
#         Contains methods
# 		to add pairs to the network from a sequence
# 		to compute scores of edges (inheritance and number of embedded trees)
#                 to reduce pairs (not used if input consists of only trees)

# '''
# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########                     AAA I/O Functions                     #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# ########
# ######## Convert Newick to a networkx Digraph with labels (and branch lengths)
# ########
# #Write length newick: convert ":" to "," and then evaluate as list of lists using ast.literal_eval
# # Then, in each list, the node is followed by the length of the incoming arc.
# # This only works as long as each branch has length and all internal nodes are labeled.
# def Newick_To_Tree(newick, current_labels = dict()):
#     newick=newick[:-1]
#     distances = False
#     #presence of : indicates the use of lengths in the trees
#     if ":" in newick:
#         distances = True
#         #taxon names may already be enclosed by " or ', otherwise, we add these now
#         if not "'" in newick and not '"' in newick:
#             newick = re.sub(r"([,\(])([a-zA-Z\d]+)", r"\1'\2", newick)
#             newick = re.sub(r"([a-zA-Z\d]):", r"\1':", newick)
#         newick = newick.replace(":",",")
#     else:
#         #taxon names may already be enclosed by " or ', otherwise, we add these now
#         if not "'" in newick and not '"' in newick:
#             newick = re.sub(r"([,\(])([a-zA-Z\d]+)", r"\1'\2", newick)
#             newick = re.sub(r"([a-zA-Z\d])([,\(\)])", r"\1'\2", newick)
#     #turn the string into a pyhton nested list using [ instead of (
#     newick = newick.replace("(","[")
#     newick = newick.replace(")","]")
#     nestedtree = ast.literal_eval(newick)
#     #parse the nested list into a list of edges with some additional information about the leaves
#     #we start with the root 2, so that we can append a root edge (1,2)
#     edges, leaves, current_labels, current_node = NestedList_To_Tree(nestedtree,2,current_labels, distances = distances)
#     #put all this information into a networkx DiGraph with or without distances/lengths
#     tree = nx.DiGraph()
#     if distances:
#         edges.append((1,2,0))
#         tree.add_weighted_edges_from(edges,weight='length')
#     else:
#         edges.append((1,2))
#         tree.add_edges_from(edges)
#     return tree, leaves, current_labels, distances

# #Auxiliary function to convert list of lists to tree (graph)
# #Works recursively, where we keep track of the nodes we have already used
# #Leaves are nodes with negative integer as ID, and already existing taxa are coupled to node IDs by current_labels.
# def NestedList_To_Tree(nestedList, next_node, current_labels, distances = False):
#     edges = []
#     leaves = set()
#     top_node = next_node
#     current_node = next_node+1
#     if distances:
#         #each element in the sublist has 2 properties, the subtree, and the length, which are adjacent in nestedList
#         for i in range(0,len(nestedList),2):
#             t = nestedList[i]
#             length = nestedList[i+1]
#             if type(t)==list: #Not a leaf
#                 edges.append((top_node,current_node,length))
#                 extra_edges, extra_leaves, current_labels, current_node = NestedList_To_Tree(t,current_node,current_labels,distances=distances)
#             else: #A leaf
#                 if str(t) not in current_labels:
#                     current_labels[str(t)] = -len(current_labels)
#                 edges.append((top_node,current_labels[str(t)],length))
#                 extra_edges = []
#                 extra_leaves = set([current_labels[str(t)]])
#             edges = edges + extra_edges
#             leaves = leaves.union(extra_leaves)
#     else:
#         #no lengths/distances, so each subtree is simply an element of nestedList
#         for t in nestedList:
#             if type(t)==list:
#                 edges.append((top_node,current_node))
#                 extra_edges, extra_leaves, current_labels, current_node = NestedList_To_Tree(t,current_node,current_labels)
#             else:
#                 if str(t) not in current_labels:
#                     current_labels[str(t)] = -len(current_labels)
#                 edges.append((top_node,current_labels[str(t)]))
#                 extra_edges = []
#                 extra_leaves = set([current_labels[str(t)]])
#             edges = edges + extra_edges
#             leaves = leaves.union(extra_leaves)
#     return edges, leaves, current_labels, current_node

# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########          AAB ANALYSIS OF SEQUENCE FOR TREES               #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# #Checks whether a given cherry-picking sequence `seq' reduces a given tree `tree'
# #if not,    returns false
# #otherwise, returns the indices of the pairs that actually reduce a cherry in the tree
# def Sequence_Reduces_Tree(seq,tree):
#     t_copy=deepcopy(tree)
#     indices = []
#     for i, pair in enumerate(seq):
#         if t_copy.reduce_pair(*pair):
#             indices +=[i]
#             if len(t_copy.nw.edges)==1:
#                 return indices
#     return False


# #Modifies a cherry-picking sequence so that it represents a network with exactly one root.
# #A sequence may be such that reconstructing a network from the sequence results in multiple roots
# #This function adds some pairs to the sequence so that the network has a single root.
# #returns the new sequence, and also modifies the sets of trees reduced by each pair in the sequence, so that the new pairs are also represented (they reduce no trees)
# def Sequence_add_roots(seq, red_trees):
#     leaves_encountered = set()
#     roots = set()
#     #The roots can be found by going back through the sequence and finding pairs where the second element has not been encountered in the sequence yet
#     for pair in reversed(seq):
#         if pair[1] not in leaves_encountered:
#             roots.add(pair[1])
#         leaves_encountered.add(pair[0])
#         leaves_encountered.add(pair[1])
#     i=0
#     roots = list(roots)
#     #Now add some pairs to make sure each second element is already part of some pair in the sequence read backwards, except for the last pair in the sequence
#     for i in range(len(roots)-1):
#        seq.append((roots[i],roots[i+1]))
#        #none of the trees are reduced by the new pairs.
#        red_trees.append(set())
#        i+=1
#     return seq, red_trees


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########         AAC INPUT SET CLASS with CPS methods              #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# #Methods for sets of phylogenetic trees
# class Input_Set:
#     def __init__(self, newick_strings = []):
#         #The dictionary of trees
#         self.trees = dict()
#         #the set of leaf labels of the trees
#         self.labels = dict()
#         self.labels_reversed = dict()
#         self.leaves = set()

#         #the current best sequence we have found for this set of trees
#         self.best_seq = None
#         #the list of reduced trees for each of the pairs in the best sequence
#         self.best_red_trees = None

#         #the best sequence for the algorithm using lengths as input as well
#         self.best_seq_with_lengths = None
#         #the sets of reduced trees for each pair in this sequence
#         self.best_seq_with_lengths_red_trees = None
#         #the height of each pair in this sequence
#         self.best_seq_with_lengths_heights = None

#         #true if distances are used
#         self.distances = True
#         #computation times
#         self.CPS_Compute_Time = 0
#         self.CPS_Compute_Reps = 0

#         #read the input trees in 'newick_strings'
#         for n in newick_strings:
#             tree = PhT()
#             self.trees[len(self.trees)] = tree
#             self.labels, distances_in_tree = tree.Tree_From_Newick(newick = n, current_labels = self.labels)
#             self.distances = self.distances and distances_in_tree
#         self.leaves = list(self.labels)

#         #make a reverse dictionary for the leaf labels, to look up the label of a given node
#         for l,i in self.labels.items():
#             self.labels_reversed[i] = l

#     #Make a deepcopy of an instance
#     def __deepcopy__(self, memodict={}):
#         copy_inputs = Input_Set()
#         copy_inputs.trees = deepcopy(self.trees,memodict)
#         copy_inputs.labels = deepcopy(self.labels,memodict)
#         copy_inputs.labels_reversed = deepcopy(self.labels_reversed,memodict)
#         copy_inputs.leaves = deepcopy(self.leaves,memodict)
# #        copy_inputs.best_seq = deepcopy(self.best_seq)
# #        copy_inputs.best_red_trees = deepcopy(self.best_red_trees)
#         return copy_inputs

#     #Find new cherry-picking sequences for the trees and update the best found
#     def CPSBound(self, repeats = 1, progress = False, track = False, lengths = False, time_limit = None):
#         #Set the specific heuristic that we use, based on the user input and whether the trees have lengths
#         Heuristic = self.CPHeuristic
#         if track and not lengths:
#             print("Tracking reducible cherries")
#             Heuristic = self.CPHeuristicStorePairs
#         if lengths:
#             if not self.distances:
#                 print("not all trees have branch lengths!")
#                 sys.exit()
#                 return []
#             print("Picking the lowest cherry")
#             Heuristic = self.CPHeuristicLengths
#             heights_best = []
#         #Initialize the recorded best sequences and corresponding data
#         best = None
#         red_trees_best = []
#         starting_time = time.time()
#         #Try as many times as required by the integer 'repeats'
#         for i in range(repeats):
#             print(i+1)
#             if lengths:
#                 new, reduced_trees, seq_heights = Heuristic(progress=progress)
#                 print("found sequence of length: "+str(len(new)))
#             else:
#                 new, reduced_trees = Heuristic(progress=progress)
#                 print("found sequence of length: "+str(len(new)))
#                 print("improving sequence")
#                 new, reduced_trees = self.Improve_Sequence(new, reduced_trees, progress = progress)
#                 print("new length = "+str(len(new)))
#             print("adding roots")
#             new,reduced_trees = Sequence_add_roots(new,reduced_trees)
#             if lengths:
#                 for i in range(len(new)-len(seq_heights)):
#                     seq_heights+=[seq_heights[-1]]
#             print("final length = "+str(len(new)))
#             if best == None or len(new) < len(best):
#                 best = new
#                 red_trees_best = reduced_trees
#                 if lengths:
#                     heights_best = seq_heights
#             print("best sequence has length "+str(len(best)))
#             self.CPS_Compute_Reps+=1
#             if time_limit and time.time()-starting_time > time_limit:
#                 break
#         self.CPS_Compute_Time += time.time() - starting_time
#         new_seq = best
#         if lengths:
#             if not self.best_seq_with_lengths or len(new_seq)<len(self.best_seq_with_lengths):
#                 converted_new_seq = []
#                 for pair in new_seq:
#                     converted_new_seq += [(self.labels_reversed[pair[0]],self.labels_reversed[pair[1]])]
#                 self.best_seq_with_lengths = converted_new_seq
#                 self.best_seq_with_lengths_red_trees = red_trees_best
#                 self.best_seq_with_lengths_heights = heights_best
#             return self.best_seq_with_lengths
#         else:
#             if not self.best_seq or len(new_seq)<len(self.best_seq):
#                 converted_new_seq = []
#                 for pair in new_seq:
#                     converted_new_seq += [(self.labels_reversed[pair[0]],self.labels_reversed[pair[1]])]
#                 self.best_seq = converted_new_seq
#                 self.best_red_trees = red_trees_best
#             return self.best_seq


#     #Version of the code that uses minimal memory: recompute reducible pairs when necessary.
#     def CPHeuristic(self, progress=False):
#         if progress:
#             print("Copying all inputs to reduce on")
#         #Works in a copy of the input trees, copy_of_inputs, because trees have to be reduced somewhere.
#         copy_of_inputs = deepcopy(self)
#         if progress:
#             print("Done, starting reduction of trees")
#         CPS = []
#         reduced_trees = []
#         candidate_leaves = deepcopy(self.leaves)
#         while copy_of_inputs.trees:
#             if progress:
#                 print("Sequence has length: "+str(len(CPS)))
#                 print(str(len(copy_of_inputs.trees))+" trees left.\n")
#                 print("Reducing trivial pairs")
#                  #First reduce trivial cherries
#             new_seq, new_red_trees = copy_of_inputs.Reduce_Trivial_Pairs(candidate_leaves)
#             if progress:
#                 print("done")
#             CPS           += new_seq
#             reduced_trees += new_red_trees
#             if len(copy_of_inputs.trees)==0:
#                 break
#              #Now reduce a random cherry from a random tree
#             random_index, random_tree = random.choice(list(copy_of_inputs.trees.items()))
#             list_of_cherries = random_tree.Find_All_Reducible_Pairs()
#             random_cherry = random.choice(list(list_of_cherries))
#             CPS         += [random_cherry]
#             reduced_by_random_cherry = copy_of_inputs.Reduce_Pair_In_All(random_cherry)
#             reduced_trees += [reduced_by_random_cherry]
#             candidate_leaves = set(random_cherry)
#         return CPS, reduced_trees

#     #Version of the code that uses more memory: stores all reducible pairs.
#     #Runs when user toggles -t or --track
#     def CPHeuristicStorePairs(self, progress=False):
#         if progress:
#             print("Copying all inputs to reduce on")
#         #Works in a copy of the input trees, copy_of_inputs, because trees have to be reduced somewhere.
#         copy_of_inputs = deepcopy(self)
#         if progress:
#             print("Done")
#         CPS = []
#         reduced_trees = []
#         candidate_leaves = deepcopy(self.leaves)
#         # Make dict of reducible pairs
#         if progress:
#             print("finding all reducible pairs")
#         reducible_pairs = self.Find_All_Pairs()
#         if progress:
#             print("found all reducible pairs")
#         while copy_of_inputs.trees:
#             if progress:
#                 print("Sequence has length: "+str(len(CPS)))
#                 print(str(len(copy_of_inputs.trees))+" trees left.\n")
#                 print("Reducing trivial pairs")
#               #First reduce trivial cherries
#             new_seq, new_red_trees, reducible_pairs = copy_of_inputs.Reduce_Trivial_Pairs_Store_Pairs(candidate_leaves, reducible_pairs)
#             if progress:
#                 print("done")
#             CPS           += new_seq
#             reduced_trees += new_red_trees
#             if len(copy_of_inputs.trees)==0:
#                 break

#                #Now reduce a random cherry from a random tree
#               #EITHER: (Get random tree, then random pair from the tree), just like in CPHeuristic
#             random_index, random_tree = random.choice(list(copy_of_inputs.trees.items()))
#             list_of_cherries = random_tree.Find_All_Reducible_Pairs()
#             random_cherry = random.choice(list(list_of_cherries))

#               #OR: (Get a random reducible pair from all pairs)
#               #Note that this would result in a different algorithm than CPHeuristic, so we use the previous option
# #            random_cherry = random.choice(list(reducible_pairs.keys()))


#             CPS         += [random_cherry]
#             #reduce all trees with this pair, this is where the list of reducible_pairs is used
#             #using the list makes it faster to find all trees that need to be reduced.
#             reduced_by_random_cherry = copy_of_inputs.Reduce_Pair_In_All(random_cherry, reducible_pairs = reducible_pairs)
#             reducible_pairs = copy_of_inputs.Update_Reducible_Pairs(reducible_pairs, reduced_by_random_cherry)
#             reduced_trees += [reduced_by_random_cherry]
#             candidate_leaves = set(random_cherry)
#         return CPS, reduced_trees


#     #Version of the code that always picks the lowest available pair
#     #Runs when user toggles -l or --lengths and all edges in the input trees have lengths.
#     def CPHeuristicLengths(self, progress=False):
#         if progress:
#             print("Copying all inputs to reduce on")
#         #Works in a copy of the input trees, copy_of_inputs, because trees have to be reduced somewhere.
#         copy_of_inputs = deepcopy(self)
#         if progress:
#             print("Done")
#         CPS = []
#         reduced_trees = []
#         heights_seq = []

#         candidate_leaves = deepcopy(self.leaves)
#         # Make dict of reducible pairs
#         if progress:
#             print("finding all reducible pairs")
#         reducible_pairs = self.Find_All_Pairs()
#         current_heights = dict()    #for each reducible pair: [0] gives height, [1] the number of trees it was computed in.

#         if progress:
#             print("found all reducible pairs")
#         while copy_of_inputs.trees:
#             if progress:
#                 print("Sequence has length: "+str(len(CPS)))
#                 print(str(len(copy_of_inputs.trees))+" trees left.\n")
#               #First reduce trivial cherries
#                 print("Reducing trivial pairs")
#             new_seq, new_red_trees, reducible_pairs, new_heights_seq = copy_of_inputs.Reduce_Trivial_Pairs_Lengths(candidate_leaves, reducible_pairs)
#             if progress:
#                 print("done")
#             CPS           += new_seq
#             reduced_trees += new_red_trees
#             heights_seq   += new_heights_seq
#             if len(copy_of_inputs.trees)==0:
#                 break

#                #Now find the lowest cherry.
#             current_heights = copy_of_inputs.Update_Heights(current_heights, reducible_pairs)
#             lowest_cherry        = None
#             lowest_height        = None
#             lowest_height_tuple  = None
#             lowest_heights_found = 1
#             for pair in reducible_pairs:
#                 height_pair_tuple = current_heights[pair][0]
#                 height_pair       = float(height_pair_tuple[0]+height_pair_tuple[1])/2
#                 new_found = False
#                 if (not lowest_height) or lowest_height>height_pair:
#                     new_found = True
#                     lowest_heights_found=1
#                 elif lowest_height==height_pair:
#                     lowest_heights_found+=1
#                     if random.random()<1/float(lowest_heights_found):
#                         new_found = True
#                 if new_found:
#                     lowest_cherry       = pair
#                     lowest_height       = height_pair
#                     lowest_height_tuple = height_pair_tuple


#             CPS         += [lowest_cherry]
#             heights_seq += [lowest_height_tuple]
#             reduced_by_lowest_cherry = copy_of_inputs.Reduce_Pair_In_All(lowest_cherry, reducible_pairs = reducible_pairs)
#             reducible_pairs = copy_of_inputs.Update_Reducible_Pairs(reducible_pairs, reduced_by_lowest_cherry)
#             reduced_trees += [reduced_by_lowest_cherry]
#             candidate_leaves = set(lowest_cherry)
#         return CPS, reduced_trees, heights_seq


#     #Returns an updated distcionary of heights of the reducible pairs
#     def Update_Heights(self,current_heights,reducible_pairs):
#         for pair, trees in reducible_pairs.items():
#             #updating is only necessary when the set of trees for that pair is changed or the reducible pair was not reducible before.
#             if not pair in current_heights or not current_heights[pair][1]==len(trees):
#                 height_pair = self.Height_Pair(pair,trees)
#                 current_heights[pair] = (height_pair, len(trees))
#         return current_heights


#     #Returns the average height of a pair in a set of trees
#     # The pair must be reducible in each tree in 'trees'
#     def Height_Pair(self, pair, trees):
#         height_pair = [0,0]
#         for t in trees:
#             height_in_t = self.trees[t].Height_Of_Cherry(*pair)
#             height_pair[0]+=height_in_t[0]
#             height_pair[1]+=height_in_t[1]
#         return [height_pair[0]/float(len(trees)),height_pair[1]/float(len(trees))]


#     #Finds the set of reducible pairs in all trees
#     #Returns a dictionary with reducible pairs as keys, and the trees they reduce as values.
#     def Find_All_Pairs(self):
#         reducible_pairs = dict()
#         for i,t in self.trees.items():
#             red_pairs_t = t.Find_All_Reducible_Pairs()
#             for pair in red_pairs_t:
#                 if pair in reducible_pairs:
#                     reducible_pairs[pair].add(i)
#                 else:
#                     reducible_pairs[pair] = set([i])
#         return reducible_pairs


#     #UReturns the updated dictionary of reducible pairs in all trees after a reduction (with the trees they reduce as values)
#     #we only need to update for the trees that got reduced: 'new_red_treed'
#     def Update_Reducible_Pairs(self, reducible_pairs, new_red_trees):
#         #Remove trees to update from all pairs
#         for pair,trees in reducible_pairs.items():
#             trees.difference_update(new_red_trees)
#             if len(trees)==0:
#                 del reducible_pairs[pair]
#         #Add the trees to the right pairs again
#         for index in new_red_trees:
#             if index in self.trees:
#                 t = self.trees[index]
#                 red_pairs_t = t.Find_All_Reducible_Pairs()
#                 for pair in red_pairs_t:
#                     if pair in reducible_pairs:
#                         reducible_pairs[pair].add(index)
#                     else:
#                         reducible_pairs[pair] = set([index])
#         return reducible_pairs

#     #reduces the given pair in all trees
#     #Returns the set of trees thet were reduced
#     #CHANGES THE SET OF TREES, ONLY PERFORM IN A COPY OF THE CLASS INSTANCE
#     def Reduce_Pair_In_All(self, pair, reducible_pairs = dict()):
#         reduced_trees_for_pair=[]
#         if pair in reducible_pairs:
#             trees_to_reduce = reducible_pairs[pair]
#         else:
#             if reducible_pairs:
#                 print("pair not found, trying all trees")
#             trees_to_reduce = deepcopy(self.trees)
#         for i in trees_to_reduce:
#             if i in self.trees:
#                 t = self.trees[i]
#                 if t.reduce_pair(*pair):
#                     reduced_trees_for_pair+=[i]
#                     if len(t.nw.edges())<=1:
#                         del self.trees[i]
#         return set(reduced_trees_for_pair)


#     #reduces the trivial pairs in the current set of trees
#     #runs efficiently by giving a set of leaves 'candidate_leaves' that may be involved in trivial pairs
#     #this set must be given; after a reduction of the pair (a,b) only using the leaves a and b works
#     #Returns the reduced pairs and the sets of trees thet were reduced
#     #CHANGES THE SET OF TREES, ONLY PERFORM IN A COPY OF THE CLASS INSTANCE
#     def Reduce_Trivial_Pairs(self,candidate_leaves):
#         seq = []
#         reduced_tree_sets = []
#         while candidate_leaves:
#             l = candidate_leaves.pop()
#             new_pairs = list(self.Trivial_Pair_With(l))
#             if new_pairs:
#                 seq += new_pairs
#                 for p in new_pairs:
#                     red_trees_p = self.Reduce_Pair_In_All(p)
#                     reduced_tree_sets += [red_trees_p]
#                     candidate_leaves   = candidate_leaves | set(p)
#         return seq, reduced_tree_sets

#     #reduces the trivial pairs in the current set of trees
#     #runs efficiently by giving a set of leaves 'candidate_leaves' that may be involved in trivial pairs
#     #this set must be given; after a reduction of the pair (a,b) only using the leaves a and b works
#     #Returns the reduced pairs and the sets of trees thet were reduced, also updates the reducible pairs.
#     #CHANGES THE SET OF TREES, ONLY PERFORM IN A COPY OF THE CLASS INSTANCE
#     def Reduce_Trivial_Pairs_Store_Pairs(self,candidate_leaves, reducible_pairs):
#         seq = []
#         reduced_tree_sets = []
#         while candidate_leaves:
#             l = candidate_leaves.pop()
#             new_pairs = list(self.Trivial_Pair_With(l))
#             if new_pairs:
# #                print("found a trivial pair")
#                 seq += new_pairs
#                 for p in new_pairs:
#                     red_trees_p = self.Reduce_Pair_In_All(p, reducible_pairs = reducible_pairs)
#                     reducible_pairs = self.Update_Reducible_Pairs(reducible_pairs, red_trees_p)
#                     reduced_tree_sets += [red_trees_p]
#                     candidate_leaves   = candidate_leaves | set(p)
#         return seq, reduced_tree_sets, reducible_pairs

#     #reduces the trivial pairs in the current set of trees with branch lengths
#     #runs efficiently by giving a set of leaves 'candidate_leaves' that may be involved in trivial pairs
#     #this set must be given; after a reduction of the pair (a,b) only using the leaves a and b works
#     #Returns the reduced pairs and the sets of trees thet were reduced, also updates the reducible pairs and their heights.
#     #CHANGES THE SET OF TREES, ONLY PERFORM IN A COPY OF THE CLASS INSTANCE
#     def Reduce_Trivial_Pairs_Lengths(self,candidate_leaves, reducible_pairs):
#         seq = []
#         reduced_tree_sets = []
#         heights_seq = []
#         while candidate_leaves:
#             l = candidate_leaves.pop()
#             new_pairs = list(self.Trivial_Pair_With(l))
#             if new_pairs:
# #                print("found a trivial pair")
#                 seq += new_pairs
#                 for p in new_pairs:
#                     height_p = self.Height_Pair(p,reducible_pairs[p])
#                     red_trees_p = self.Reduce_Pair_In_All(p, reducible_pairs = reducible_pairs)
#                     heights_seq += [height_p]
#                     reducible_pairs = self.Update_Reducible_Pairs(reducible_pairs, red_trees_p)
#                     reduced_tree_sets += [red_trees_p]
#                     candidate_leaves   = candidate_leaves | set(p)
#         return seq, reduced_tree_sets, reducible_pairs, heights_seq

#     #Returns all trivial pairs involving the leaf l
#     def Trivial_Pair_With(self,l):
#         pairs = set()
#         #Go through all trees t with index i.
#         for i,t in self.trees.items():
#             #If the leaf occurs in t
#             if l in t.leaves:
#                 #Compute reducible pairs of t with the leaf as first coordinate
#                 pairs_in_t = t.Find_Pairs_With_First(l)
#                 #If we did not have a set of candidate pairs yet, use pairs_in_t
#                 if not pairs:
#                     pairs = pairs_in_t
#                 #Else, the candidate pairs must also be in t, so take intersection
#                 else:
#                     pairs = pairs & pairs_in_t
#                 #If we do not have any candidate pairs after checking a tree with l as leaf, we stop.
#                 if not pairs:
#                     break
#         return pairs


#     #Improves a sequence 'CPS' for the input trees by removing elements and checking whether the new sequence still reduces all trees
#     #Returns this improved sequence and the corresponding sets of reduced trees for each pair.
#     def Improve_Sequence(self, CPS, reduced_trees, progress = False):
#         seq = deepcopy(CPS)
#         i=0
#         while i<len(seq):
#             redundant = True
#             relevant_tree_indices = reduced_trees[i]
#             new_relevant_pairs_for_trees = dict()
#             for j in relevant_tree_indices:
#                 #Check if the shorter sequence reduces the trees, and if so, record which pairs reduced a cherry in which tree
#                 new_relevant_pairs_for_trees[j] = Sequence_Reduces_Tree(seq[:i]+seq[i+1:],self.trees[j])
#                 if not new_relevant_pairs_for_trees[j]:
#                     redundant = False
#                     break
#             if redundant:
#                 #Remove the ith element from seq and reduced trees
#                 seq.pop(i)
#                 reduced_trees.pop(i)
#                 #Update reduced_trees for the relevant trees
#                    #First remove all j in relevant_tree_indices
#                 for tree_index_set in reduced_trees:
#                     for j in relevant_tree_indices:
#                         tree_index_set.discard(j)
#                    #Now add them back at the right places, according to "new_relevant_pairs_for_trees"
#                 for j in relevant_tree_indices:
#                     for index in new_relevant_pairs_for_trees[j]:
#                         reduced_trees[index].add(j)
#                 if progress:
#                     print("New length is "+str(len(seq)))
#                     print("Continue at position "+str(i))
#             else:
#                 i+=1
#         return seq, reduced_trees


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########              AAD PHYLOGENETIC TREE CLASS                  #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# #A class representing a phylogenetic tree
# #Contains methods to reduce trees
# class PhT:
#     def __init__(self):
#         #the actual graph
#         self.nw = nx.DiGraph()
#         #the set of leaf labels of the network
#         self.leaves = set()

#     #Builds a tree from a newick string
#     def Tree_From_Newick(self, newick = None, current_labels = dict()):
#         self.nw, self.leaves, current_labels, distances = Newick_To_Tree(newick, current_labels)
#         return current_labels, distances

#     #Checks whether the pair (x,y) forms a cherry in the tree
#     def Is_Cherry(self,x,y):
#         if (not x in self.leaves) or (not y in self.leaves):
#             return False
#         px=-1
#         py=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         for p in self.nw.predecessors(y):
#             py=p
#         return px==py

#     #Returns the height of (x,y) if it is a cherry:
#     #     i.e.: length(p,x)+length(p,y)/2
#     #Returns false otherwise
#     def Height_Of_Cherry(self,x,y):
#         if (not x in self.leaves) or (not y in self.leaves):
#             return False
#         px=-1
#         py=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         for p in self.nw.predecessors(y):
#             py=p
#         if px==py:
#             height = [float(self.nw[px][x]['length']),float(self.nw[py][y]['length'])]
#             return height
#         return False


#     #suppresses a degree-2 node v and returns true if successful
#     # the new arc has length length(p,v)+length(v,c)
#     #returns false if v is not a degree-2 node
#     def Clean_Node(self,v):
#         if self.nw.out_degree(v)==1 and self.nw.in_degree(v)==1:
#             pv=-1
#             for p in self.nw.predecessors(v):
#                 pv=p
#             cv=-1
#             for c in self.nw.successors(v):
#                 cv=c
#             self.nw.add_edges_from([(pv,cv,self.nw[pv][v])])
#             if 'length' in self.nw[pv][v] and 'length' in self.nw[v][cv]:
#                 self.nw[pv][cv]['length']=self.nw[pv][v]['length']+self.nw[v][cv]['length']
#             self.nw.remove_node(v)
#             return True
#         return False

#     #reduces the pair (x,y) in the tree if it is present as cherry
#     # i.e., removes the leaf x and its incoming arc, and then cleans up its parent node.
#     # note that if px, and py have different lengths, the length of px is lost in the new network.
#     #returns true if successful and false otherwise
#     def reduce_pair(self,x,y):
#         if not x in self.leaves or not y in self.leaves:
#             return False
#         px=-1
#         py=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         for p in self.nw.predecessors(y):
#             py=p
#         if self.Is_Cherry(x,y):
#             self.nw.remove_node(x)
#             self.leaves.remove(x)
#             self.Clean_Node(py)
#             return True
#         return False


#     #Returns all reducible pairs in the tree involving x, where x is the first element
#     def Find_Pairs_With_First(self,x):
#         pairs = set()
#         px=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         if self.nw.out_degree(px)>1:
#             for cpx in self.nw.successors(px):
#                 if cpx in self.leaves:
#                     pairs.add((x,cpx))
#         return pairs-set([(x,x)])

#     #Returns all reducible pairs in the tree
#     def Find_All_Reducible_Pairs(self):
#         red_pairs=set()
#         for l in self.leaves:
#             red_pairs=red_pairs.union(self.Find_Pairs_With_First(l))
#         return red_pairs


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########                     AAE CutTree CLASS                     #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


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


# ################################################################################
# ################################################################################
# ################################################################################
# ########                                                           #############
# ########              AAF PHYLOGENETIC NETWORK CLASS               #############
# ########                                                           #############
# ################################################################################
# ################################################################################
# ################################################################################


# #A class for phylogenetic networks
# class PhN:
#     def __init__(self, seq = None, newick = None, best_tree_from_network = None, reduced_trees = None, heights = None):
#         #the actual graph
#         self.nw = nx.DiGraph()
#         #the set of leaf labels of the network
#         self.leaves = set()
#         #a dictionary giving the node for a given leaf label
#         self.labels = dict()
#         #the number of nodes in the graph
#         self.no_nodes = 0
#         self.leaf_nodes = dict()
#         self.TCS=seq
#         self.CPS=seq
#         self.newick=newick
#         self.reducible_pairs=set()
#         self.reticulated_cherries=set()
#         self.cherries=set()
#         self.level = None
#         self.no_embedded_trees = 0
#         #if a cherry-picking sequence is given, build the network from this sequence
#         if seq:
#             total_len = len(seq)
#             current_trees_embedded = set()
#             #Creates a phylogenetic network from a cherry picking sequence:
#             if reduced_trees:
#                 for i,pair in enumerate(reversed(seq)):
#                     if heights:
#                         self.add_pair(*pair,red_trees = reduced_trees[total_len-1-i], current_trees=current_trees_embedded, height = heights[total_len-1-i])
#                         current_trees_embedded = current_trees_embedded | reduced_trees[total_len-1-i]
#                     else:
#                         self.add_pair(*pair,red_trees = reduced_trees[total_len-1-i], current_trees=current_trees_embedded)
#                 self.no_embedded_trees = len(current_trees_embedded)
#                 self.ScoreEdges()
#             else:
#                 for pair in reversed(seq):
#                     self.add_pair(*pair)
#         #if a newick string is given, build the network from the newick string
#         elif newick:
#             self.newick = newick
#             network, self.leaves, self.labels = Newick_To_Network(newick)
#             self.nw = network
#             self.no_nodes = len(list(self.nw))
#             self.Compute_Leaf_Nodes()
#         #if a network 'best_tree_from_network' is given, extract the best tree from this network and use this tree as the network
#         elif best_tree_from_network:
#             self.nw.add_edges_from(best_tree_from_network.Best_Tree())
#             self.labels = best_tree_from_network.labels
#             self.leaf_nodes = best_tree_from_network.leaf_nodes
#             self.leaves = best_tree_from_network.leaves
#             self.no_nodes = best_tree_from_network.no_nodes
#             #self.Clean_Up()


#     #Checks if the graph in the instance is a network
#     # i.e., checks whether it has a single root, all levaes are labeled, there are no degree-2 nodes, and each non-leaf non-root node is either a reticulation or a tree node.
#     def IsANetwork(self):
#         rootFound = False
#         for v in self.nw.nodes:
#             if self.nw.in_degree(v) == 0:
#                 if rootFound:
#                     print("Multiple Roots!")
#                     return False
#                 rootFound = True
#             if self.nw.out_degree(v) == 0 and not v in self.leaf_nodes:
#                 print("Unlabelled leaf!")
#                 return False
#             if self.nw.out_degree(v) == 1 and self.nw.in_degree(v) == 1:
#                 print("Suppressable node, do not forget to clean up your network!")
#                 return False
#             if self.nw.out_degree(v) > 1 and self.nw.in_degree(v) > 1:
#                 print("Combined split and reticulation!")
#                 return False
#         return True


#     #Returns the leaf nodes of the network
#     def Compute_Leaf_Nodes(self):
#         self.leaf_nodes = dict()
#         for v in self.labels:
#             self.leaf_nodes[self.labels[v]]=v

#     #Returns the level (retirculations in the largest biconnected component) of the network
#     def Level(self, recompute = False):
#         if recompute or not self.level:
#             blobs = nx.biconnected_component_edges(self.nw.to_undirected())
#             lvl = 0
#             for b_edges in blobs:
#                 blob = nx.Graph(b_edges)
#                 lvl = max(lvl, blob.number_of_edges() - blob.number_of_nodes() +1)
#             self.level = lvl
#         return self.level


#     #Add the probabilities for each edge, returns true
#     # probabilties of hybrid edges are stored as edge property 'prob' (i.e. fraction of TREES THAT GO THROUGH THE RETICULATION NODE that use this edge)
#     # probabilities of all edges are stored in 'probability_all'
#     # the fraction of ALL INPUT TREES going through the edge is stored in 'frac_of_trees'
#     def ScoreEdges(self):
#         for node in self.nw.nodes:
#             #for the hybrid edges
#             if self.nw.in_degree(node)>1:
#                 total_trees = 0.0
#                 for parent in self.nw.predecessors(node):
#                    total_trees+= self.nw[parent][node]['no_of_trees']
#                 for parent in self.nw.predecessors(node):
#                    self.nw[parent][node]['prob'] = self.nw[parent][node]['no_of_trees']/total_trees
#                    self.nw[parent][node]['probability_all'] = self.nw[parent][node]['prob']
#                    if self.no_embedded_trees > 0:
#                        self.nw[parent][node]['frac_of_trees'] = self.nw[parent][node]['no_of_trees']/float(self.no_embedded_trees)
#             #and also for the non-hybrid edges (their probability is one)
#             else:
#                 for parent in self.nw.predecessors(node):
#                     self.nw[parent][node]['probability_all'] = 1
#                     if self.no_embedded_trees > 0:
#                         self.nw[parent][node]['frac_of_trees'] = self.nw[parent][node]['no_of_trees']/float(self.no_embedded_trees)
#         return True


#     #Adds a pair to the network, using the construction from a cherry-picking sequence
#     # returns false if y is not yet in the network and the network is not empty
#     def add_pair(self,x,y, red_trees = set(), current_trees = set(), height = [1,1]):
#         #if the network is empty, create a cherry (x,y)
#         if len(self.leaves)==0:
#             self.nw.add_edge(0,1,no_of_trees=len(red_trees), length=0)
#             self.nw.add_edge(1,2,no_of_trees=len(red_trees), length=height[0])
#             self.nw.add_edge(1,3,no_of_trees=len(red_trees), length=height[1])
#             self.leaves = set([x,y])
#             self.labels[x]=2
#             self.labels[y]=3
#             self.leaf_nodes[2]=x
#             self.leaf_nodes[3]=y
#             self.no_nodes=4
#             return True
#         #if y is not in the network return false, as there is no way to add the pair and get a phylogenetic network
#         if y not in self.leaves:
#             return False
#         #add the pair to the existing network
#         node_y=self.labels[y]
#         parent_node_y = -1
#         for p in self.nw.predecessors(node_y):
#             parent_node_y=p

#         #first add all edges around y
#         length_incoming_y = self.nw[parent_node_y][node_y]['length']
#         no_of_trees_incoming_y = self.nw[parent_node_y][node_y]['no_of_trees']
#         height_goal_x = height[0]
#         if height[1]<length_incoming_y:
#             height_pair_y_real = height[1]
#         else:
#             height_pair_y_real = length_incoming_y
#             height_goal_x+=height[1]-height_pair_y_real


#         self.nw.add_edge(node_y,self.no_nodes,no_of_trees=no_of_trees_incoming_y+len(red_trees-current_trees), length=height_pair_y_real)
#         self.nw[parent_node_y][node_y]['length'] = length_incoming_y - height_pair_y_real
#         self.leaf_nodes.pop(self.labels[y],False)
#         self.labels[y]=self.no_nodes
#         self.leaf_nodes[self.no_nodes]=y

#         #Now also add edges around x
#         #x is not yet in the network, so make a cherry (x,y)
#         if x not in self.leaves:
#             self.nw.add_edge(node_y,self.no_nodes+1,no_of_trees=len(red_trees), length=height_goal_x)
#             self.leaves.add(x)
#             self.labels[x]=self.no_nodes+1
#             self.leaf_nodes[self.no_nodes+1]=x
#             self.no_nodes+=2
#         #x is already in the network, so create a reticulate cherry (x,y)
#         else:
#             node_x=self.labels[x]
#             for parent in self.nw.predecessors(node_x):
#                 px = parent
#             length_incoming_x = self.nw[px][node_x]['length']
#             no_of_trees_incoming_x = self.nw[px][node_x]['no_of_trees']
#             #if x is below a reticulation, and the height of the new pair is above the height of this reticulation, add the new hybrid arc to the existing reticulation
#             if self.nw.in_degree(px)>1 and length_incoming_x<=height_goal_x:
#                 self.nw.add_edge(node_y,px,no_of_trees=len(red_trees), length=height_goal_x-length_incoming_x)
#                 self.nw[px][node_x]['no_of_trees']+=len(red_trees)
#                 self.no_nodes+=1
#             #create a new reticulation vertex above x to attach the hybrid arc to
#             else:
#                 height_pair_x = min(height_goal_x,length_incoming_x)
#                 self.nw.add_edge(node_y,node_x,no_of_trees=len(red_trees), length=height_goal_x-height_pair_x)
#                 self.nw.add_edge(node_x,self.no_nodes+1,no_of_trees = no_of_trees_incoming_x+len(red_trees), length = height_pair_x)
#                 self.nw[px][node_x]['length'] = length_incoming_x - height_pair_x
#                 self.leaf_nodes.pop(self.labels[x],False)
#                 self.labels[x]=self.no_nodes+1
#                 self.leaf_nodes[self.no_nodes+1]=x
#                 self.no_nodes+=2
#         return True


#     #returns the 'best tree' in the network
#     # i.e., pick the best incoming arc for each reticulation, and use this to find a tree.
#     def Best_Tree(self):
#         edges = []
#         for v in self.nw.nodes:
#             best_p=None
#             best_value=-1
#             for p in self.nw.predecessors(v):
#                 if self.nw[p][v]['no_of_trees']>best_value:
#                     best_value = self.nw[p][v]['no_of_trees']
#                     best_p=p
#             if best_p:
#                 edges.append((best_p,v,self.nw[best_p][v]))
#         return edges

#     #Returns a subnetwork by selecting the highest scoring hybrid arcs (we pick 'reticulations' of these) and extending this to a network, so that the resulting network
#     # the threshold may be the fraction of trees that uses the edge, or the probability of the edge.
#     def SelectSubNetworkByReticulations(self,type_is_probability = True, reticulations = 0):
#         restrictedNetwork = PhN(best_tree_from_network = self)
#         score_type = 'frac_of_trees'
#         if type_is_probability:
#             score_type = 'probability_all'
#         for i in range(reticulations):
#             best_score = -1
#             best_edge = None
#             for e in self.nw.edges:
#                 if (not restrictedNetwork.nw.has_edge(e[0],e[1])) and self.nw[e[0]][e[1]][score_type]>best_score:
#                     best_edge = e
#                     best_score = self.nw[e[0]][e[1]][score_type]
#             #If we have already selected the whole network
#             if not best_edge:
#                 break
#             restrictedNetwork.nw.add_edges_from([(best_edge[0],best_edge[1],self.nw[best_edge[0]][best_edge[1]])])
#             up_node = best_edge[1]
#             while restrictedNetwork.nw.in_degree(up_node) == 0:
#                 best_score = -1
#                 best_parent = None
#                 for parent in self.nw.predecesors(up_node):
#                     if self.nw[parent][up_node][score_type] > best_score:
#                         best_score = self.nw[parent][up_node][score_type]
#                         best_parent = parent
#                 restrictedNetwork.nw.add_edges_from([(best_parent,up_node,self.nw[best_parent][up_node])])
#                 up_node = best_parent
#         restrictedNetwork.Clean_Up()
#         #Debug: check if we find a valid network
#         restrictedNetwork.IsANetwork()
#         restrictedNetwork.ScoreEdges()
#         return restrictedNetwork

#     #Returns a subnetwork by selecting all hybrid arcs that have a score above a given threshold and extending this to a network
#     # the threshold may be the fraction of trees that uses the edge, or the probability of the edge.
#     def SelectSubNetworkByScore(self,type_is_probability = True, score = .5):
#         restrictedNetwork = PhN(best_tree_from_network = self)
#         score_type = 'frac_of_trees'
#         if type_is_probability:
#             score_type = 'probability_all'
#         done = False
#         while not done:
#             best_score = -1
#             best_edge = None
#             for e in self.nw.edges:
#                 if (not restrictedNetwork.nw.has_edge(e[0],e[1])) and self.nw[e[0]][e[1]][score_type]>best_score:
#                     best_edge = e
#                     best_score = self.nw[e[0]][e[1]][score_type]
#             if best_score >= score:

#                 restrictedNetwork.nw.add_edges_from([(best_edge[0],best_edge[1],self.nw[best_edge[0]][best_edge[1]])])
#                 up_node = best_edge[1]
#                 while restrictedNetwork.nw.in_degree(up_node) == 0:
#                     best_score = -1
#                     best_parent = None
#                     for parent in self.nw.predecesors(up_node):
#                         if self.nw[parent][up_node][score_type] > best_score:
#                             best_score = self.nw[parent][up_node][score_type]
#                             best_parent = parent
#                     restrictedNetwork.nw.add_edges_from([(best_parent,up_node,self.nw[best_parent][up_node])])
#                     up_node = best_parent
#             else:
#                 done = True
#         restrictedNetwork.Clean_Up()
#         #Debug: check if we find a valid network
#         restrictedNetwork.IsANetwork()
#         restrictedNetwork.ScoreEdges()
#         return restrictedNetwork


#     #Returns true if (x_label,y_label) forms a cherry in the network, false otherwise
#     def Is_Cherry(self,x_label,y_label):
#         if not x_label in self.leaves or not x_label in self.leaves:
#             return False
#         x = self.labels[x_label]
#         y = self.labels[y_label]
#         px=-1
#         py=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         for p in self.nw.predecessors(y):
#             py=p
#         return px==py

#     #Returns true if (x_label,y_label) forms a reticulate cherry in the network, false otherwise
#     def Is_Ret_Cherry(self,x_label,y_label):
#         if not x_label in self.leaves or not x_label in self.leaves:
#             return False
#         x = self.labels[x_label]
#         y = self.labels[y_label]
#         px=-1
#         py=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         for p in self.nw.predecessors(y):
#             py=p
#         return (self.nw.in_degree(px)>1) and self.nw.out_degree(px)==1 and (py in self.nw.predecessors(px))

#     #Returns true if (x_label,y_label) forms a reducible pair in the network, false otherwise
#     def Is_Reducible_Pair(self,x_label,y_label):
#         if not x_label in self.leaves or not x_label in self.leaves:
#             return False
#         x = self.labels[x_label]
#         y = self.labels[y_label]
#         px=-1
#         py=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         for p in self.nw.predecessors(y):
#             py=p
#         return (px==py) or (self.nw.in_degree(px)>1 and self.nw.out_degree(px)==1 and py in self.nw.predecessors(px))

#     #Iteratively removes unlabeled leaves until none are left, then supresses al degree-2 nodes
#     def Clean_Up(self):
#         nodes_to_check = set(deepcopy(self.nw.nodes))
#         while nodes_to_check:
#             v = nodes_to_check.pop()
#             if self.nw.out_degree(v)==0 and v not in self.leaf_nodes:
#                 for p in self.nw.predecessors(v):
#                     nodes_to_check.add(p)
#                 self.nw.remove_node(v)
#         list_nodes = deepcopy(self.nw.nodes)
#         for v in list_nodes:
#             self.Clean_Node(v)

#     #supresses v if it is a degree-2 node
#     def Clean_Node(self,v):
#         if self.nw.out_degree(v)==1 and self.nw.in_degree(v)==1:
#             pv=-1
#             for p in self.nw.predecessors(v):
#                 pv=p
#             cv=-1
#             for c in self.nw.successors(v):
#                 cv=c
#             self.nw.add_edges_from([(pv,cv,self.nw[v][cv])])
#             if self.nw[pv][v]['length'] and self.nw[v][cv]['length']:
#                 self.nw[pv][cv]['length'] = self.nw[pv][v]['length'] + self.nw[v][cv]['length']
#             self.nw.remove_node(v)
#             return True
#         return False

#     #reduces the pair (x_label,y_label) if it is reducible in the network
#     #returns a new set reducible pairs that involve the leaves x_label and y_label
#     def reduce_pair(self,x_label,y_label):
#         if not x_label in self.leaves or not y_label in self.leaves:
#             return set()
#         x = self.labels[x_label]
#         y = self.labels[y_label]
#         px=-1
#         py=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         for p in self.nw.predecessors(y):
#             py=p
#         if self.Is_Cherry(x_label,y_label):
#             self.reducible_pairs.difference_update(set([(x_label,y_label),(y_label,x_label)]))
#             self.nw.remove_node(x)
#             self.leaves.remove(x_label)
#             self.labels.pop(x_label,False)
#             self.Clean_Node(py)
#             #AddCherriesInvolving y
#             new_pairs = set([("no_leaf","no_leaf")])|self.Find_Pairs_With_First(y_label)| self.Find_Pairs_With_Second(y_label)
#             self.reducible_pairs=self.reducible_pairs.union(new_pairs-set([("no_leaf","no_leaf")]))
#             return new_pairs
#         if self.Is_Ret_Cherry(x_label,y_label):
#             self.reducible_pairs.difference_update(set([(x_label,y_label),(y_label,x_label)]))
#             self.nw.remove_edge(py,px)
#             self.Clean_Node(px)
#             self.Clean_Node(py)
#             #AddCherriesInvolving x and y
#             new_pairs = set([("no_leaf","no_leaf")])|self.Find_Pairs_With_First(x_label) |self.Find_Pairs_With_Second(x_label)|self.Find_Pairs_With_First(y_label)|self.Find_Pairs_With_Second(y_label)
#             self.reducible_pairs=self.reducible_pairs.union(new_pairs-set([("no_leaf","no_leaf")]))
#             return new_pairs
#         return set()


#     #Returns all reducible pairs in the network where x_label is the first element of the pair
#     def Find_Pairs_With_First(self,x_label):
#         pairs = set()
#         x = self.labels[x_label]
#         px=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         if self.nw.in_degree(px)>1:
#             for ppx in self.nw.predecessors(px):
#                 for cppx in self.nw.successors(ppx):
#                     if cppx in self.leaf_nodes:
#                         pairs.add((x_label,self.leaf_nodes[cppx]))
#         if self.nw.out_degree(px)>1:
#             for cpx in self.nw.successors(px):
#                 if cpx in self.leaf_nodes:
#                     pairs.add((x_label,self.leaf_nodes[cpx]))
#         return pairs-set([(x_label,x_label)])

#     #Returns all reducible pairs in the network where x_label is the second element of the pair
#     def Find_Pairs_With_Second(self,x_label):
#         pairs = set()
#         x = self.labels[x_label]
#         px=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         if self.nw.out_degree(px)>1:
#             for cpx in self.nw.successors(px):
#                 if cpx in self.leaf_nodes:
#                     pairs.add((self.leaf_nodes[cpx],x_label))
#                 if self.nw.in_degree(cpx)>1:
#                     for ccpx in self.nw.successors(cpx):
#                         if ccpx in self.leaf_nodes:
#                             pairs.add((self.leaf_nodes[ccpx],x_label))
#         return pairs-set([(x_label,x_label)])

#     #Returns all cherries and reticulate cherries in the network where x_label is the second element of the pair
#     def Find_Pairs_With_Second_Separated(self,x_label):
#         cherries = set()
#         retic_cherries = set()
#         x = self.labels[x_label]
#         px=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         if self.nw.out_degree(px)>1:
#             for cpx in self.nw.successors(px):
#                 if cpx in self.leaf_nodes:
#                     cherries.add((self.leaf_nodes[cpx],x_label))
#                 if self.nw.in_degree(cpx)>1:
#                     for ccpx in self.nw.successors(cpx):
#                         if ccpx in self.leaf_nodes:
#                             retic_cherries.add((self.leaf_nodes[ccpx],x_label))
#         return cherries-set([(x_label,x_label)]), retic_cherries-set([(x_label,x_label)])

#     #Returns all reticulate cherries in the network where x_label is the second element of the pair
#     def Find_Retic_Cherry_Second(self,x_label):
#         pairs = set()
#         x = self.labels[x_label]
#         px=-1
#         for p in self.nw.predecessors(x):
#             px=p
#         if self.nw.out_degree(px)>1:
#             for cpx in self.nw.successors(px):
#                 if self.nw.in_degree(cpx)>1:
#                     for ccpx in self.nw.successors(cpx):
#                         if ccpx in self.leaf_nodes:
#                             pairs.add((self.leaf_nodes[ccpx],x_label))
#         return pairs-set([(x_label,x_label)])

#     #Returns all reducible pairs in the network
#     # also sets this as value of self.reducible_pairs
#     def Find_All_Reducible_Pairs(self):
#         self.reducible_pairs=set()
#         for l in self.leaves:
#             self.reducible_pairs=self.reducible_pairs.union(self.Find_Pairs_With_First(l))
#             self.reducible_pairs=self.reducible_pairs.union(self.Find_Pairs_With_Second(l))
#         return self.reducible_pairs

#     #Returns all cherries and reticulate cherries in the network
#     # also sets these as values of self.cherries and self.reticulate_cherries
#     def Find_All_Reducible_Pairs_Separated(self):
#         self.reticulated_cherries=set()
#         self.cherries=set()
#         for l in self.leaves:
#             added_cherries, added_retic_cherries = self.Find_Pairs_With_Second_Separated(l)
#             self.cherries=self.reticulated_cherries.union(added_cherries)
#             self.reticulated_cherries=self.reticulated_cherries.union(added_retic_cherries)
#         return(self.cherries, self.reticulated_cherries)
