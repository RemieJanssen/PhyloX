import os
import sys
import math
import re
import ast
import random
import numpy as np
import pandas as pd
from pathlib import Path
import multiprocessing
from multiprocessing import Manager
from ExtractNetworkProperties import *
import itertools




def Newick_To_Network(newick):
    newick=newick[:-1]

    #remove internal labels
    newick = re.sub(r"I([\d]+)", "", newick)

    #remove lengths
    newick = re.sub(r"([\d]+)\.([\d]+)", "", newick)
    newick = re.sub(r"E-[\d]+", "", newick)
    newick = re.sub(r":", "", newick)

    #make into list format
    newick = newick.replace("(","[")
    newick = newick.replace(")","]")
    newick = re.sub(r"\]\#H([\d]+)", r",#R\1]", newick)
    newick = re.sub(r"#([RH])([\d]+)", r"'#\1\2'", newick)


    #add "" if necessary    
    newick = re.sub(r"([ABCD])", r"'\1'", newick)
    newick = re.sub(r" ", "", newick)
    
    
    nestedtree = ast.literal_eval(newick)
    edges, leaves, label_set, current_node = NestedList_To_Tree(nestedtree,1)
    edges.append([0,1])
    ret_labels = dict()
    leaf_labels = dict()
    for l in leaves:
        if len(l)>2 and (l[:2]=="#H" or l[:2]=="#R"):
            ret_labels[l[2:]]=[]
        else:
            leaf_labels[l]=[]
    for l in label_set:
        if len(l[0])>2 and (l[0][:2]=="#H" or l[0][:2]=="#R"):
            if l[0][1]=='H':
                ret_labels[l[0][2:]]+=[l[1]]
            else:
                ret_labels[l[0][2:]]=[l[1]]+ret_labels[l[0][2:]]
        else:
            leaf_labels[l[0]]+=[l[1]]
    network = nx.DiGraph()        
    network.add_edges_from(edges)
    for retic in ret_labels:
        r = ret_labels[retic]
        receiving = r[0]
        parent_receiving = 0
        for p in network.predecessors(receiving):
            parent_receiving = p
        network.remove_node(receiving)
        for v in r[1:]:
            network.add_edge(v,parent_receiving)
            network = nx.contracted_edge(network,(v,parent_receiving))
            network.remove_edge(v,v)
            parent_receiving = v
    leaves = set()
    for l in leaf_labels:
         leaf_labels[l]=leaf_labels[l][0]
         leaves.add(l)
    return network, leaves, leaf_labels
    
    
def NestedList_To_Tree(nestedList,next_node):
    edges = []
    leaves = set()
    labels = []
    top_node = next_node
    current_node = next_node+1
    for t in nestedList:
        edges.append((top_node,current_node))
        if type(t)==list: 
            extra_edges, extra_leaves, extra_labels, current_node = NestedList_To_Tree(t,current_node)
        else: 
            extra_edges = []
            extra_leaves = set([str(t)])
            extra_labels = [[str(t), current_node]]
            current_node+=1
        edges = edges + extra_edges
        leaves = leaves.union(extra_leaves)
        labels = labels + extra_labels
    return edges, leaves, labels, current_node




network,_,leaf_labels = Newick_To_Network(newick_string)
reverse_labels = {x:y for (y,x) in leaf_labels.items()}


