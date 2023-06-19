import os
import sys
import math
import random
import numpy as np
import networkx as nx
from pathlib import Path
import multiprocessing
from multiprocessing import Manager
from BetaSplittingNetwork_Tools import *


############################################
# MAIN
############################################

############################################
# General parameters.

# Desired tree size.
leaves          = [8,16,32]
# Range of reticulations.
reticulations   = [1,2,3,4,5,6,7,8,9,10,15,20,25,30,40,50,100,200,300,400,500,1000,2000,5000]
# Edge adding methods.
methods         = ["Horizontal"]#[None,0.01,0.02,0.10,0.20]
# Number of networks per tree, reticulation number, and method.
nb_networks     = 100

# Folder 
folder = "balanced_and_caterpillar/"

############################################


def Balanced_Tree(n):
    tree = nx.DiGraph()
    tree.add_nodes_from(range(1,n+1))
    roots = list(range(1,n+1))
    curr = n+1
    while len(roots)>1:
        while len(roots)>1:
            r1 = roots.pop()
            r2 = roots.pop()
            tree.add_edges_from([(curr,r1),(curr,r2)])
            curr+=1
        roots = [x for x in tree.nodes() if tree.in_degree(x)==0]
    return tree



def Caterpillar_Tree(n):
    tree = nx.DiGraph()
    tree.add_node(n+1)
    for i in range(1,n):
        tree.add_edges_from([(n+i,i),(n+i,n+i+1)])
    tree=nx.relabel_nodes(tree,{2*n:n})
    return tree             

# Set all parameters for network simulation.
all_sim = []
for retics in reticulations:
    for method in methods:
        for index in range(nb_networks):
            all_sim+=[(retics,method,index)]

for out_type in ["el","pl"]:
    this_path = folder+out_type
    if not os.path.exists(this_path):
        os.makedirs(this_path) 
    

# Simulate networks
for tree_type in ["balanced","caterpillar"]:
    for n in leaves:
        print(tree_type,n)
        if tree_type =="balanced":
            tree = Balanced_Tree(n)
        else:
            tree = Caterpillar_Tree(n)
        tree_prefix = "/"+tree_type+"_"+str(n)
        for out_type in ["el","pl"]:
            this_path = folder+out_type+tree_prefix+"/0"
            if not os.path.exists(this_path):
                os.makedirs(this_path) 
            f=open(this_path+"/1.txt", 'w')
            f.write(OutputNetwork(tree,out_type))
            f.close()    

        
        # Method to simulate a single network.
        def runsubproc(pars):
            # Simulate tree.
            retics, method, index  = pars
            network  = GenerateNetwork(tree,retics,method)
            filepath_network  = "/"+str(retics)+"/"+str(method)
            index_txt      = "/"+str(index+1)+".txt"
            for out_type in ["el", "pl"]:
                this_path = folder+out_type+tree_prefix+filepath_network
                if not os.path.exists(this_path):
                    os.makedirs(this_path) 
                filename = this_path+index_txt
                if not os.path.exists(filename):
                    f=open(filename, 'w')
                    f.write(OutputNetwork(network,out_type))
                    f.close()    

        #p = multiprocessing.Pool(nb_threads)
        #p.map(runsubproc, all_sim)
        for par in all_sim:
            runsubproc(par)
        


