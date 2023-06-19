# import CPH
# import os
# import sys
# import csv

# errorInFile = False

# #####################################################################
# #####################################################################
# ##############          Terminal arguments                ###########
# #####################################################################
# #####################################################################

# #Set default values
# option_out = False
# option_out_argument = ""
# option_file = False
# option_file_argument = ""
# option_progress = False
# option_bestTree = False
# option_repeats = False
# option_track = False
# option_repeats_argument = 1
# option_lengths = False
# option_timeLimit = False
# option_timeLimit_argument = None
# option_scoreSubnetwork = False
# option_scoreSubnetwork_argument = None
# option_reticsSubnetwork = False
# option_reticsSubnetwork_argument = None
# option_help = False

# #Read the arguments
# i = 1
# while i < len(sys.argv):
#     arg= sys.argv[i]
#     if arg == "-f" or arg == "--file":
#         option_file = True
#         i+=1
#         option_file_argument = sys.argv[i]
#         option_out_argument = option_file_argument+"OUT"
#     if arg == "-o" or arg == "--output":
#         option_out = True
#         i+=1
#         option_out_argument = sys.argv[i]
#     if arg == "-p" or arg == "--progress":
#         option_progress = True
#     if arg == "-t" or arg == "--track":
#         option_track = True
#     if arg == "-l" or arg == "--lengths":
#         option_lengths = True
#     if arg == "-r" or arg == "--repeats":
#         option_repeats = True
#         i+=1
#         option_repeats_argument = int(sys.argv[i])
#     if arg == "-rsn" or arg == "--reticssubnetwork":
#         option_reticsSubnetwork = True
#         i+=1
#         option_reticsSubnetwork_argument = int(sys.argv[i])
#     if arg == "-ssn" or arg == "--scoresubnetwork":
#         option_scoreSubnetwork = True
#         i+=1
#         option_scoreSubnetwork_argument = float(sys.argv[i])
#     if arg == "-tl" or arg == "--timelimit":
#         option_timeLimit = True
#         i+=1
#         option_timeLimit_argument = float(sys.argv[i])
#     if arg == "-bt" or arg == "--besttree":
#         option_bestTree = True
#     if arg == "-h" or arg == "--help":
#         option_help = True
#     i += 1

# #Output the help text to the terminal if no argument is given, or if the help option is chosen.
# if len(sys.argv)==1 or option_help:
#     print("Mandatory arguments:\n -f or --file followed by the input file in simple Newick format. \n\nOptional arguments:\n -o or --output followed by the output file containing the sequence, Default: [input_filename]+OUT\n -p or --progress for extra progress reports during the running of the algorithm\n -r or --repeats followed by an integer, the number of repeated attempts to find a good sequence.\n -tl or --timelimit followed by a float number of seconds, the timelimit given to the algorithm. The code will finish a repetition it is currently executing.\n -bt or --besttree for extracting the best tree out of the network, i.e. selecting the most used edges per reticulation node.\n -t or --track For tracking all reducible pairs during the algorithm.\n -l or --lengths For picking the lowest reducible cherry every time. If there are multiple, pick one at random.\n -rsn or --reticssubnetwork followed by an integer, the reticulation number of the subnetwork, to extract a subnetwork with this number of reticulations by iteratively adding reticulations with high probability (if probabilities are given) or by adding the reticulation edge that is used most by all input trees. \n -ssn or --scoresubnetwork followed by a threshold number (float) to extract a subnetwork by using all edges with a score above this threshold. If the network has probabilities, this score is the probability of an edge, otherwise it is the fraction of the input trees that uses this edge in their embedding.")
#     sys.exit()


# #####################################################################
# #####################################################################
# ##############             Read the input                 ###########
# #####################################################################
# #####################################################################

# #Empty set of inputs
# inputs = []

# #Read each line of the input file with name set by "option_file_argument"
# f = open("./"+option_file_argument, "rt")
# reader = csv.reader(f, delimiter='~', quotechar='|')
# for row in reader:
#     inputs.append(str(row[0]))
# f.close()

# #Make the set of inputs usable for all algorithms: use the CPH class
# tree_set = CPH.Input_Set(newick_strings = inputs)


# #####################################################################
# #####################################################################
# ##############            Find the sequence               ###########
# #####################################################################
# #####################################################################

# #Run the heuristic to find a cherry-pciking sequence `seq' for the set of input trees.
# #Arguments are set as given by the terminal arguments
# seq = tree_set.CPSBound(repeats = option_repeats_argument,
#                     progress    = option_progress,
#                     track       = option_track,
#                     lengths     = option_lengths,
#                     time_limit  = option_timeLimit_argument)


# #####################################################################
# #####################################################################
# ##############                   Output                   ###########
# #####################################################################
# #####################################################################


# #Open the output file for writing
# f= open(option_out_argument,"w+")

# #Output the raw network.
# f.write("Best network in Newick format;\r\n")
# if option_lengths:
#     network = CPH.PhN(seq = tree_set.best_seq_with_lengths,
#             reduced_trees = tree_set.best_seq_with_lengths_red_trees,
#             heights       = tree_set.best_seq_with_lengths_heights)
# else:
#     network = CPH.PhN(seq = tree_set.best_seq,
#             reduced_trees = tree_set.best_red_trees)
# #Cut the network into a tree to find the newick string
# cuttree = CPH.CutTree(network    = network.nw,
#                     current_node = network.no_nodes,
#                     leaf_labels  = network.leaf_nodes)
# f.write(cuttree.Newick(probabilities = True)+"\r\n")
# f.write(cuttree.Newick(probabilities = False)+"\r\n")

# #Output the desired subnetworks.
# f.write("Good subnetworks in Newick format;\r\n")
# #Find the best tree inside the network
# if option_bestTree:
#     resNW = network.SelectSubNetworkByReticulations(type_is_probability=True, reticulations = 0)
#     resNW.Clean_Up()
#     cuttree = CPH.CutTree(network = resNW.nw, current_node = resNW.no_nodes, leaf_labels  = resNW.leaf_nodes)
#     f.write(cuttree.Newick(probabilities = True)+"\r\n")
#     f.write(cuttree.Newick(probabilities = False)+"\r\n")
# #Find the best subnetwork with the given number of reticulations: option_reticsSubnetwork_argument from option -rsn
# if option_reticsSubnetwork:
#     #probabilities on edges
#     resNW = network.SelectSubNetworkByReticulations(type_is_probability=True, reticulations = option_reticsSubnetwork_argument)
#     cuttree = CPH.CutTree(network = resNW.nw, current_node = resNW.no_nodes, leaf_labels = resNW.leaf_nodes)
#     f.write(cuttree.Newick(probabilities = True)+"\r\n")
#     f.write(cuttree.Newick(probabilities = False)+"\r\n")
#     #no probabilities on edges
#     resNW = network.SelectSubNetworkByReticulations(type_is_probability=False, reticulations = option_reticsSubnetwork_argument)
#     cuttree = CPH.CutTree(network = resNW.nw, current_node = resNW.no_nodes, leaf_labels = resNW.leaf_nodes)
#     f.write(cuttree.Newick(probabilities = True)+"\r\n")
#     f.write(cuttree.Newick(probabilities = False)+"\r\n")
# #Find the best subnetwork with the given number of reticulations: option_reticsSubnetwork_argument from option -rsn
# if option_scoreSubnetwork:
#     #probabilities on edges
#     resNW = network.SelectSubNetworkByScore(type_is_probability=True, score = option_scoreSubnetwork_argument)
#     cuttree = CPH.CutTree(network = resNW.nw, current_node = resNW.no_nodes, leaf_labels = resNW.leaf_nodes)
#     f.write(cuttree.Newick(probabilities = True)+"\r\n")
#     f.write(cuttree.Newick(probabilities = False)+"\r\n")
#     #no probabilities on edges
#     resNW = network.SelectSubNetworkByScore(type_is_probability=False, score = option_scoreSubnetwork_argument)
#     cuttree = CPH.CutTree(network = resNW.nw, current_node = resNW.no_nodes, leaf_labels = resNW.leaf_nodes)
#     f.write(cuttree.Newick(probabilities = True)+"\r\n")
#     f.write(cuttree.Newick(probabilities = False)+"\r\n")


# #Output additional information about the given arguments
# f.write("END\r\n")
# f.write("Selected by score: \r\n")
# f.write(str(option_scoreSubnetwork_argument)+"\r\n")
# f.write("Selected by reticulation number: \r\n")
# f.write(str(option_reticsSubnetwork_argument)+"\r\n\r\n")

# #Output the computation time for the heuristic
# f.write("Computation time; CPS Heuristic reps\r\n")
# f.write(str(tree_set.CPS_Compute_Time) +";" + str(tree_set.CPS_Compute_Reps)+"\r\n\r\n")

# #Output the raw sequence found by the heuristic
# last = '0'
# f.write("Best sequence\r\n")
# f.write("Pair; reduced trees; (height of reduced pair)\r\n")
# for i, pair in enumerate(seq):
#     if option_lengths:
#         f.write("("+pair[0]+","+pair[1]+"); "+str(tree_set.best_seq_with_lengths_red_trees[i])+"); "+str(tree_set.best_seq_with_lengths_heights[i])+"\r\n")
#     else:
#         f.write("("+pair[0]+","+pair[1]+"); "+str(tree_set.best_red_trees[i])+"\r\n")
#     last = pair[1]
# f.write("("+last+",-)\r\n")

# #close the output file
# f.close()
