############################################
#
# Adapted from BetasplittingModel.py
# which contains the code for Beta-splitting model (Aldous1996)
# Written for python3
#
############################################

import os
import sys
import math
import random
import numpy as np
import networkx as nx
from scipy.special import gamma, loggamma



############################################
# Simulation functions
############################################

# a_n is a normalizing constant defined in 
# Equation (2) of Aldous1996 (so the sum of 
# the values is equal to 1. It is not 
# computed here to save time.
def a_n(beta):
    return 1

# Compute the "probability" to split n in (i|n-1), where i=1,..,n-1
def computeProb(n,beta):
    q_n = []
    for i in range(1,n):
        q_i_n = np.exp((loggamma(beta+i+1)+loggamma(beta+n-i+1))-((a_n(beta)+loggamma(i+1)+loggamma(n-i+1))))
        q_n.append(q_i_n)
    return q_n

# n: number of tips
def simulateBetaSplitting(n, beta):
    # Initialize tree.
    tree = nx.DiGraph()
    tree.add_node(n+1)
    tree.node[n+1]['label'] = n
    last_internal_node       = n+1
    last_leaf_node           = 0 
    queue                    = [n+1]
    # Insert one node at each iteration.
    while queue:
        node        = queue.pop()
        n_node      = tree.node[node].get('label')
        # Internal node. Splits again.
        if (n_node > 1): 
            # Compute the "probability" to split n in (i|n-1), where i=1,..,n-1
            q_n   = computeProb(n_node,beta)
            split = random.choices(population=list(range(1,n_node)), weights=q_n, k=1)[0]
            # Create children.
            for new_n in [split,n_node-split]:
                if new_n == 1:
                    tree.add_edge(node,last_leaf_node+1)
                    last_leaf_node+=1
                else:
                    tree.add_edge(node,last_internal_node+1)
                    tree.node[last_internal_node+1]['label']=new_n
                    queue.append(last_internal_node+1)
                    last_internal_node+=1
    # Return tree.
    return tree


def GenerateNetwork(tree,r,method):
    network = tree.copy()
    if method=="Horizontal":
        leaves = [x for x in tree.nodes() if tree.out_degree(x)==0]
    for i in range(r):
        if method==None:
            AddEdgeUniform(network)
        elif method=="Horizontal":
            AddEdgeHorizontal(network,leaves=leaves)
        else:
            AddEdgeLocal(network,stop_prob=method)
    return network


def AddEdgeBetween(network,edge1,edge2,new_nodes=None):
    #if we dont have new nodes yet, determine new nodes
    if new_nodes==None:
        max_node = max(network.nodes())
        new_nodes = (max_node+1,max_node+2)
    #make sure edge2 is not above edge1
    if nx.has_path(network,edge2[1],edge1[0]):
        edge1,edge2=edge2,edge1
    length1 = network[edge1[0]][edge1[1]].get('length')
    prob1 = network[edge1[0]][edge1[1]].get('prob')    
    length2 = network[edge2[0]][edge2[1]].get('length')
    prob2 = network[edge2[0]][edge2[1]].get('prob')    
    #add an edge from edge1 to edge2
    network.remove_edges_from([edge1,edge2])
    network.add_edges_from([(edge1[0],new_nodes[0]),(new_nodes[0],edge1[1]),(edge2[0],new_nodes[1]),(new_nodes[1],edge2[1]),(new_nodes[0],new_nodes[1])])
    
    #If the original network had lengths and probabilities, add these to the subdivided arcs as well.
    if length1!=None:
        l11 = random.random()*length1
        l12 = length1-l11
        network[edge1[0]][new_nodes[0]]['length']=l11
        network[new_nodes[0]][edge1[1]]['length']=l12
    if prob1!=None:
        network[new_nodes[0]][edge1[1]]['prob']=prob1
    if length2!=None:
        l21 = random.random()*length2
        l22 = length2-l21
        network[edge2[0]][new_nodes[1]]['length']=l21
        network[new_nodes[1]][edge2[1]]['length']=l22
    if prob2!=None:
        network[new_nodes[1]][edge2[1]]['prob']=prob2
    #TODO: add probabilities and length to the other new arcs as well?

    

#Pick two edges uniformly at random and add an edge between these
def AddEdgeHorizontal(network,leaves = None,new_nodes=None):
    if leaves==None:
        leaves = [x for x in network.nodes() if network.out_degree(x)==0]
    leaf_indices = np.random.choice(range(len(leaves)),2,replace=False)
    l0 = leaves[leaf_indices[0]]
    l1 = leaves[leaf_indices[1]]
    e0 = list(network.in_edges(l0))[0]
    e1 = list(network.in_edges(l1))[0]
    AddEdgeBetween(network,e0,e1,new_nodes=new_nodes)


#Pick two edges uniformly at random and add an edge between these
def AddEdgeUniform(network,new_nodes=None):
    edges = list(network.edges())
    edge_indices = np.random.choice(range(len(edges)),2,replace=False)
    AddEdgeBetween(network,edges[edge_indices[0]],edges[edge_indices[1]],new_nodes=new_nodes)


    
    
#Pick one edge, move a random number of edges through the network to find a second edge
#Add and edge between the two edges.    
def AddEdgeLocal(network,new_nodes=None,stop_prob=0.2,max_steps=None,max_tries=None):
    try_number = 1
    while max_tries==None or try_number<=max_tries:
        #Pick a random edge
        edge1 = random.choice(list(network.edges()))
        edge2 = None
        #Initiate the random walk, by choosing an orientation
        previous_node = random.choice(edge1)
        current_node  = edge1[0]
        if current_node == previous_node:
            current_node = edge1[1]
        #Take a number of steps
        step_number =1
        while max_steps==None or step_number<=max_steps:
            previous_node,current_node = current_node,random.choice(list(network.successors(current_node))+list(network.predecessors(current_node)))
            if random.random()<stop_prob:
                break
            step_number+=1
        #Set the new edge
        edge2 = (previous_node,current_node)
        if edge2 not in network.edges():
            edge2 = (current_node,previous_node)
        #Add an edge if possible, otherwise repeat the search
        if edge1!=edge2:
            break
        try_number+=1
    AddEdgeBetween(network,edge1,edge2,new_nodes=new_nodes)


def OutputNetwork(network,out_type):
    output = ""
    if out_type=="el":
        for i,e in enumerate(network.edges):
            if i!=0:
                output+="\r\n"
            output += str(e[0])+" "+str(e[1])
    elif out_type=="pl":
        for i,v in enumerate(network.nodes):
            if i!=0:
                output+="\r\n"
            output += str(v)
            for p in network.predecessors(v):
                output+=" "+str(p)
    return output



