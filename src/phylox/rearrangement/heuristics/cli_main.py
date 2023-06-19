# # These functions are implementations of the algorithms of Remie Janssen's PhD thesis

# from __future__ import print_function
# import networkx as nx
# from RearrDistance_Tools import *
# import ast
# import os
# import sys
# import re
# import time

# ##PARAMETERS
# filename = None
# edges = False
# tailMoves = True
# headMoves = True
# bfSearch = False
# simple=False
# time_limit=False
# repeats=1
# runtime=False
# deep_dive=False
# randomNodes=False
# ##


# ###############################2. I/O############################
# option_help = False
# i = 1
# while i < len(sys.argv):
#     arg= sys.argv[i]
#     if arg == "-f" or arg == "--filename":
#         i+=1
#         filename = str(sys.argv[i])
#     if arg == "-tl" or arg == "--timelimit":
#         i+=1
#         time_limit = float(sys.argv[i])
#     if arg == "-e" or arg == "--edges":
#         edges = True
#     if arg == "-t" or arg == "--tail":
#         headMoves = False
#     if arg == "-h" or arg == "--head":
#         tailMoves = False
#     if arg == "-bfs" or arg == "--bfsearch":
#         bfSearch = True
#     if arg == "-dd" or arg == "--deepdive":
#         deep_dive = True
#     if arg == "-s" or arg == "--simple":
#         simple = True
#     if arg == "-rt" or arg == "--runtime":
#         runtime = True
#     if arg == "-r" or arg == "--random":
#         randomNodes = True
#     if arg == "-rep" or arg == "--repeats":
#         i+=1
#         repeats = int(sys.argv[i])
#     if arg == "-help" or arg == "--help":
#         option_help = True
#     i += 1

# if len(sys.argv)==1 or option_help or (deep_dive and bfSearch):
#     print("Mandatory arguments:\n -f or --filename followed by the filename of the file containing two networks \n\nOptional arguments:\n -e or --edges if the input file contains a list of edges in the form [(x1,y1),...,(xn,yn)] with xi and yi integers or strings in the form \"string\". If this option is not selected, the input is assumed to consist of two newick strings.\n -t or --tail for only using tail moves, instead of tail and head moves.\n -h or --head for only using head moves, instead of tail and head moves.\n -r or --random to pick random nodes whenever arbitrary nodes are required in the heuristics.\n -rt or --runtime to also record the runtime of the algorithm and the time to read the trees from file.\n -bfs or --bfsearch for using a breadth first search to find the an optimal sequence.\n -tl or --timelimit followed by a number of seconds to set a timelimit for the bfs. If no answer is found before the time runs out, the algorithm just returns a lower bound on the distance.\n\nThe output is given as a list of moves in the format:\n  moving_edge, moving_endpoint, to_edge")
#     sys.exit()


# ####################################################
# ####################################################
# ####################################################
# #############                          #############
# #############           MAIN           #############
# #############                          #############
# ####################################################
# ####################################################
# ####################################################


# if runtime:
#     time_start_reading = time.time()

# test = open(filename, "r")
# line1 = test.read()
# line1 = line1.split("\n")
# test.close()
# if edges:
#     N = nx.DiGraph()
#     M = nx.DiGraph()
#     N.add_edges_from(ast.literal_eval(line1[0]))
#     M.add_edges_from(ast.literal_eval(line1[1]))
#     rootN=Root(N)
#     if N.out_degree(rootN)==2:
#        N.add_edges_from([('rho',rootN)])
#     rootM=Root(M)
#     if M.out_degree(rootM)==2:
#        M.add_edges_from([('rho',rootM)])
#     N = NetworkLeafToLabel(N)
#     M = NetworkLeafToLabel(M)
#     label_attribute=None
# else:
#     N = Newick_To_Network(line1[0])
#     M = Newick_To_Network(line1[1])
#     if not simple:
#         print("The networks as list of edges, with node names as used in the computed sequence of moves")
#         print("Network 1:")
#         for e in N.edges():
#             label = N.node[e[1]].get('label')
#             if label:
#                 print(str(e[0])+" "+str(e[1])+" = leaf: "+str(label))
#             else:
#                 print(str(e[0])+" "+str(e[1]))
#         print("Network 2:")
#         for e in M.edges():
#             label = M.node[e[1]].get('label')
#             if label:
#                 print(str(e[0])+" "+str(e[1])+" = leaf: "+str(label))
#             else:
#                 print(str(e[0])+" "+str(e[1]))


# if runtime:
#     time_start_searching = time.time()


# if bfSearch:
#     if not simple:
#         print("Computing a shortest sequence using breadth first search.")
#         #Note, the code uses a DFS with incremented depth bound
#         #Otherwise, the queue gets too large for memory
#     sequence = Depth_First(N,M,tail_moves=tailMoves,head_moves=headMoves,max_time=time_limit,show_bounds=(not simple))
# elif deep_dive:
#     if not simple:
#         print("Computing a shortest sequence by passing the depth first search tree once, choosing the best branch (according to the GL algorithm) each time.")
#         #Note, the code uses a DFS with incremented depth bound
#         #Otherwise, the queue gets too large for memory
#     sequence = Deep_Dive_Scored(N,M,head_moves=headMoves)


# else:
#     if not tailMoves:
#         if not simple:
#             print("Computing a sequence using the `red-line' heuristic.")
#         if randomNodes:
#             if not simple:
#                 print("Picking random nodes whenever arbitrary nodes are required.")
#             sequence = Red_Line_Random(N,M,repeats=repeats)
#         else:
#             sequence = Red_Line(N,M)
#     else:
#         if not simple:
#             print("Computing a sequence using the `green-line' heuristic.")
#         if randomNodes:
#             if not simple:
#                 print("Picking random nodes whenever arbitrary nodes are required.")
#             sequence = Green_Line_Random(N,M,head_moves=headMoves,repeats=repeats)
#         else:
#             sequence = Green_Line(N,M,head_moves=headMoves)

# if runtime:
#     time_stop_searching = time.time()
#     if simple:
#         print(";"+str(time_start_searching-time_start_reading)+";"+str(time_stop_searching-time_start_searching),end='')
#     else:
#         print("Reading the networks cost "+str(time_start_searching-time_start_reading)+" seconds.")
#         print("Finding the sequence cost "+str(time_stop_searching-time_start_searching)+" seconds.")


# if not simple:
#     print("Sequence:")
# if sequence==False:
#     if not simple:
#         print("There is no sequence between the networks.")
#     else:
#         print(";None;None",end='')
#     sys.exit()
# if type(sequence)==int:
#     if simple:
#         print(";>="+str(sequence)+";?",end='')
#     else:
#         print("No network was found within the time limit.")
#         print("The distance between the networks is at least: "+str(sequence))
#     sys.exit()
# if edges:
#     sequence = ReplaceNodeNamesByOriginal(N,sequence)

# #Print the output


# if simple:
#     print(";"+str(len(sequence))+";"+str(sequence),end='')
# else:
#     for move in sequence:
#         if len(move)==4:
#             print(str(move[0])+" "+str(move[1])+" "+str(move[3]))
#         if len(move)==3:
#             print(str(move[0])+" "+str(move[1])+" "+str(move[2]))
