import networkx as nx
from RearrDistance_Tools import *
from BetaSplittingNetwork_Tools import *
import os
import sys


"""
folderName = "2Leaf012Retic"
leaves  = 2
max_retics = 2
index_offset  = 0
beta    = -1.0
type_probs = {"horizontal":0.8,"vertPlus":0.1,"vertMin":0.1}
correct_symmetries = True

#period between samples
burn_in             = 500
#number of sampled networks
number_of_networks  = 5000
"""



folderName = "3Leaf012Retic"
leaves  = 3
max_retics = 2
index_offset  = 0
beta    = -1.0
type_probs = {"horizontal":0.8,"vertPlus":0.1,"vertMin":0.1}
correct_symmetries = True

#period between samples
burn_in             = 1000
#number of sampled networks
number_of_networks  = 10000




i=0
while i < len(sys.argv):
    arg= sys.argv[i]
    if arg == "-nw" or arg == "--basenetwork":
        i+=1
        arg1 = str(sys.argv[i])
        if arg1[0]=="-":
            filename = arg1
        else:
            leaves      = int(sys.argv[i])
            beta        = float(sys.argv[i+1])
            max_retics  = int(sys.argv[i+2])
    if arg == "-n" or arg == "--networks":
        i+=1
        number_of_networks = int(sys.argv[i])
    if arg == "-b" or arg == "--burnin":
        i+=1
        burn_in = int(sys.argv[i])
    i += 1


for out_type in ["el","pl"]:
    this_path = "./"+folderName+"_"+str(beta)+"/"+out_type+"/"
    if not os.path.exists(this_path):
        os.makedirs(this_path) 



def AutomorphismCount(network):#TODO
    print("not implemented yet")
    return 1
    



def AcceptanceProb(move,leaves,curr_retics,max_retics,type_probs,symmetries=False):
    p=0
    if move[0]=="horizontal":
        p = 1
    if move[0]=="vertPlus":        
        no_edges_network = float(2*leaves + 3*curr_retics - 1)
        no_edges_network_after  = no_edges_network+3
        p = (type_probs["vertMin"]/type_probs["vertPlus"])  *  no_edges_network**2/(no_edges_network_after)
#        p = (type_probs["vertMin"]/type_probs["vertPlus"])*(((no_edges_network)**2 * (max_retics-curr_retics)**2)/(no_edges_network+3))
#        if symmetries:
#            p*= 1/float(max_retics-curr_retics)**2
    if move[0]=="vertMin":
        no_edges_network        = float(2*leaves + 3*curr_retics - 1)
        no_edges_network_after  = no_edges_network-3
        if no_edges_network>3:
            p = (type_probs["vertPlus"]/type_probs["vertMin"]) * no_edges_network/(no_edges_network_after**2)
#            p = (type_probs["vertPlus"]/type_probs["vertMin"])*(no_edges_network/((no_edges_network_after)**2 * (max_retics-(curr_retics-1))**2))
#            if symmetries:
#                p*= float(max_retics-(curr_retics-1))**2
    if symmetries:
        #TODO: correct for number of representations, i.e., symmetries.
        p*=1
    return p



network = simulateBetaSplitting(leaves, beta)
curr_retics = 0        
root = Root(network)
network.add_edges_from([(0,root)])
print(network.nodes())
available_reticulations     = ["r"+str(k+1) for k in range(max_retics)]
available_tree_nodes        = [str(k+2*leaves) for k in range(max_retics)]

for index in range(number_of_networks):
    non_moves = 0
    start_time = time.time()
    for j in range(burn_in):
        move = RandomMove(network,available_tree_nodes=available_tree_nodes,available_reticulations=available_reticulations,type_probs=type_probs)
        result = False
        if random.random() < AcceptanceProb(move,leaves,curr_retics,max_retics,type_probs,symmetries=correct_symmetries):
            if not (move[0]=="vertPlus" and curr_retics==max_retics):
                result = DoMove(network,move)                
        if result:
            if move[0]=="horizontal":
                network = result
            if move[0]=="vertPlus":
                curr_retics +=1
                available_tree_nodes.remove(move[3])
                available_reticulations.remove(move[4])
                network = result
            if move[0]=="vertMin":
                curr_retics -=1
                available_tree_nodes+=[move[1][0]]
                available_reticulations+=[move[1][1]]
                network = result
        else:
            non_moves+=1
    print(curr_retics, index+index_offset, non_moves, time.time()-start_time)
    for out_type in ["el","pl"]:
        this_path = "./"+folderName+"_"+str(beta)+"/"+out_type+"/"+str(max_retics)+"_"+str(index+index_offset)+".txt"
        f=open(this_path, 'w')
        f.write(OutputNetwork(network,out_type))
        f.close()

        

