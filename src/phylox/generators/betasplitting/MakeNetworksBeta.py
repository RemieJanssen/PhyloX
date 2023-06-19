import os
import sys
import math
import random
import numpy as np
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
leaves          = 10
# Number of simulated trees per tree size.
nb_trees        = 50
# Range of beta parameter.
betas           = [0.0,-0.5,-1.0,-1.5]
# Range of reticulations.
reticulations   = [6,7,8,9,10,11,12,13,14,15,16,17,18,19]#[1,2,3,4,5,10,20,50,100,200,500]
# Edge adding methods.
methods         = ["Horizontal",None,0.002,0.020,0.200]
# Number of networks per tree, reticulation number, and method.
nb_networks     = 10
# Number of cores in your computer (trees are simulated in parallel to save time).
nb_threads      = 4

# Folder 
folder = "Beta_Splitting_Networks_"+str(leaves)+"/"

############################################



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
for b in betas:
    for t in range(nb_trees):
        print(b,t+1)
        #Check whether the tree exists:
        tree_prefix  = "/"+str(b)+"/"+str(t+1)
        tree_path = folder+"el"+tree_prefix+"/0/1.txt"
        if os.path.exists(tree_path):
            tree        = nx.read_edgelist(tree_path,create_using=nx.DiGraph()) 
            relabel_dict = {node: int(node) for node in tree.nodes}
            tree=nx.relabel_nodes(tree,relabel_dict)             
        else:
            tree        = simulateBetaSplitting(leaves, b)
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
            network                = GenerateNetwork(tree,retics,method)
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
        


