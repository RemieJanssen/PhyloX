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
from BetaSplittingNetwork_Tools import *



############################################
# I/O
############################################
    # General parameters and default values
# All desired tree sizes.
tree_size       = 100
# Range of beta parameter.
beta            = -0.0
# Number of reticulations
retics          = 10
# Edge addition type; 
#  default = None, in which case edges are added uniformly; 
#  otherwise, float value gives the probability of stopping the random walk
local           = None
# output type: pl=parent list, el=edge list, nw=newick (not implemented)
out_type        = "el"


option_help = False
i = 1
while i < len(sys.argv):
    arg= sys.argv[i]
    if arg == "-n" or arg == "--number_of_tips":
        i+=1
        tree_size = int(sys.argv[i])
    if arg == "-b" or arg == "--beta":
        i+=1
        beta = float(sys.argv[i])
    if arg == "-r" or arg == "--reticulations":
        i+=1
        retics = int(sys.argv[i])
    if arg == "-l" or arg == "--local":
        i+=1
        local = float(sys.argv[i])
    if arg == "-pl" or arg == "--parent_list":
        out_type = "pl"
    if arg == "-el" or arg == "--edge_list":
        out_type = "el"
    if arg == "-h" or arg == "--help":
        option_help = True
    i += 1


if option_help:
    print("Generates a tree under the beta-splitting model (Aldous 1996) and then adds edges to generate a network. The edges are added one by one, by randomly choosing two edges in the network and attaching a new edge between them. The two edges can be randomly chosen in two ways: uniformly at random [Default]; or by chosing one edge uniformly at random and then performing a random walk to find a second edge. The second option has a parameter `stop_prob' that determines the length of the random walk: after each step of the walk, the walk stops with probability `stop_prob. Hence, the higher this parameter, the more local the extra reticulation edges become.\n\nOptional arguments:\n -n or --number_of_tips followed by an int [Default=100] to set the number of tips of the network.\n -b or --beta followed by a float [Default=0.0] to set the beta parameter of the tree generating beta-splitting model.\n -r or --reticulations followed by an int [Default=10] to set the number of reticulations (i.e., extra edges).\n -l or --local followed by a float to use the local edge adding method. The float is the stop probability for the random walk.\n -pl or --parent_list to get output as a parent list.\n -el or --edge_list to get output as an edge list.")
    sys.exit()




############################################
# MAIN: Simulate network
############################################



tree = simulateBetaSplitting(tree_size, beta)
network = GenerateNetwork(tree,retics,local)
print(OutputNetwork(network,out_type))

