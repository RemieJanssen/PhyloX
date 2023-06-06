import networkx as nx
from RearrDistance_Tools import *
from BetaSplittingNetwork_Tools import *
import os
import sys


"""
leaves  = 50
retics  = [1,2,3,4,5,10,20,50,100,200,500]
folderName = "n_"+str(leaves)+"_k_"+"_".join(map(str,retics))
index_offset  = 0
betas   = [-1.0]#[-1.0,-1.5]
type_probs = {"horizontal":1.0,"vertPlus":0.0,"vertMin":0.0}
"""


leaves  = 10
retics  = list(range(21))
folderName = "n_"+str(leaves)+"_k_"+"_".join(map(str,retics))
index_offset  = 0
betas   = [-1.0]#[-1.0,-1.5]
type_probs = {"horizontal":1.0,"vertPlus":0.0,"vertMin":0.0}


#period between samples
burn_in             = 500
#number of sampled networks
number_of_networks  = 500


i=0
while i < len(sys.argv):
    arg= sys.argv[i]
    if arg == "-nw" or arg == "--basenetwork":
        i+=1
        arg1 = str(sys.argv[i])
        if arg1[0]=="-":
            filename = arg1
        else:
            leaves = int(sys.argv[i])
            beta   = float(sys.argv[i+1])
            retics = [int(sys.argv[i+2])]
    if arg == "-n" or arg == "--networks":
        i+=1
        number_of_Networks = int(sys.argv[i])
    if arg == "-b" or arg == "--burnin":
        i+=1
        burn_in = int(sys.argv[i])
    i += 1


for beta in betas:
    for out_type in ["el","pl"]:
        this_path = "./"+folderName+"_"+str(beta)+"/"+out_type+"/"
        if not os.path.exists(this_path):
            os.makedirs(this_path) 

"""
def RandomMove(network):
    edges = list(network.edges())
    edge_indices = np.random.choice(range(len(edges)),2,replace=False)
    moving_edge = edges[edge_indices[0]]
    to_edge = edges[edge_indices[1]]
    endpoint = random.choice(moving_edge)
    return (moving_edge,endpoint,to_edge)
"""

for beta in betas:
    for r in retics:
        tree = simulateBetaSplitting(leaves, beta)
        network = GenerateNetwork(tree,r,None)
        root = Root(network)
        network.add_edges_from([(max(map(int,network.nodes()))+1,root)])
        for index in range(number_of_networks):
            non_moves = 0
            start_time = time.time()
            for j in range(burn_in):
    #            candidate_moves = AllValidMoves(network,tail_moves=True,head_moves=True)
    #            move = random.choice(candidate_moves)

                move = RandomMove(network,type_probs = type_probs)
                result = DoMove(network,move)
                if result:
                    network = result
                else:
                    non_moves+=1
            print(r, index+index_offset, non_moves, time.time()-start_time)
            for out_type in ["el","pl"]:
                this_path = "./"+folderName+"_"+str(beta)+"/"+out_type+"/"+str(r)+"_"+str(index+index_offset)+".txt"
                f=open(this_path, 'w')
                f.write(OutputNetwork(network,out_type))
                f.close()

#    print(OutputNetwork(network,"el"))
#    print(*[(int(x[0]),int(x[1])) for x in list(network.edges())])
        

