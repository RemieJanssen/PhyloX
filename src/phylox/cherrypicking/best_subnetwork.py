# # Add the probabilities for each edge, returns true
# # probabilties of hybrid edges are stored as edge property 'prob' (i.e. fraction of TREES THAT GO THROUGH THE RETICULATION NODE that use this edge)
# # probabilities of all edges are stored in 'probability_all'
# # the fraction of ALL INPUT TREES going through the edge is stored in 'frac_of_trees'
# def ScoreEdges(network):
#     for node in network.nodes:
#         # for the hybrid edges
#         if network.in_degree(node) > 1:
#             total_trees = 0.0
#             for parent in network.predecessors(node):
#                 total_trees += network[parent][node]["no_of_trees"]
#             for parent in network.predecessors(node):
#                 network[parent][node]["prob"] = (
#                     network[parent][node]["no_of_trees"] / total_trees
#                 )
#                 network[parent][node]["probability_all"] = network[parent][node]["prob"]
#                 if network.no_embedded_trees > 0:
#                     network[parent][node]["frac_of_trees"] = network[parent][node][
#                         "no_of_trees"
#                     ] / float(network.no_embedded_trees)
#         # and also for the non-hybrid edges (their probability is one)
#         else:
#             for parent in network.predecessors(node):
#                 network[parent][node]["probability_all"] = 1
#                 if network.no_embedded_trees > 0:
#                     network[parent][node]["frac_of_trees"] = network[parent][node][
#                         "no_of_trees"
#                     ] / float(network.no_embedded_trees)
#     return True


# # returns the 'best tree' in the network
# # i.e., pick the best incoming arc for each reticulation, and use this to find a tree.
# def Best_Tree(network):
#     edges = []
#     for v in network.nodes:
#         best_p = None
#         best_value = -1
#         for p in network.predecessors(v):
#             if network[p][v]["no_of_trees"] > best_value:
#                 best_value = network[p][v]["no_of_trees"]
#                 best_p = p
#         if best_p:
#             edges.append((best_p, v, network[best_p][v]))
#     return edges


# # Returns a subnetwork by selecting the highest scoring hybrid arcs (we pick 'reticulations' of these) and extending this to a network, so that the resulting network
# # the threshold may be the fraction of trees that uses the edge, or the probability of the edge.
# def SelectSubNetworkByReticulations(network, type_is_probability=True, reticulations=0):
#     restrictedNetwork = PhN(best_tree_from_network=network)
#     score_type = "frac_of_trees"
#     if type_is_probability:
#         score_type = "probability_all"
#     for i in range(reticulations):
#         best_score = -1
#         best_edge = None
#         for e in network.edges:
#             if (not restrictedNetwork.has_edge(e[0], e[1])) and network[e[0]][e[1]][
#                 score_type
#             ] > best_score:
#                 best_edge = e
#                 best_score = network[e[0]][e[1]][score_type]
#         # If we have already selected the whole network
#         if not best_edge:
#             break
#         restrictedNetwork.add_edges_from(
#             [(best_edge[0], best_edge[1], network[best_edge[0]][best_edge[1]])]
#         )
#         up_node = best_edge[1]
#         while restrictedNetwork.in_degree(up_node) == 0:
#             best_score = -1
#             best_parent = None
#             for parent in network.predecesors(up_node):
#                 if network[parent][up_node][score_type] > best_score:
#                     best_score = network[parent][up_node][score_type]
#                     best_parent = parent
#             restrictedNetwork.add_edges_from(
#                 [(best_parent, up_node, network[best_parent][up_node])]
#             )
#             up_node = best_parent
#     restrictedNetwork.Clean_Up()
#     # Debug: check if we find a valid network
#     restrictedNetwork.IsANetwork()
#     restrictedNetwork.ScoreEdges()
#     return restrictedNetwork


# # Returns a subnetwork by selecting all hybrid arcs that have a score above a given threshold and extending this to a network
# # the threshold may be the fraction of trees that uses the edge, or the probability of the edge.
# def SelectSubNetworkByScore(network, type_is_probability=True, score=0.5):
#     restrictedNetwork = PhN(best_tree_from_network=network)
#     score_type = "frac_of_trees"
#     if type_is_probability:
#         score_type = "probability_all"
#     done = False
#     while not done:
#         best_score = -1
#         best_edge = None
#         for e in network.edges:
#             if (not restrictedNetwork.has_edge(e[0], e[1])) and network[e[0]][e[1]][
#                 score_type
#             ] > best_score:
#                 best_edge = e
#                 best_score = network[e[0]][e[1]][score_type]
#         if best_score >= score:

#             restrictedNetwork.add_edges_from(
#                 [(best_edge[0], best_edge[1], network[best_edge[0]][best_edge[1]])]
#             )
#             up_node = best_edge[1]
#             while restrictedNetwork.in_degree(up_node) == 0:
#                 best_score = -1
#                 best_parent = None
#                 for parent in network.predecesors(up_node):
#                     if network[parent][up_node][score_type] > best_score:
#                         best_score = network[parent][up_node][score_type]
#                         best_parent = parent
#                 restrictedNetwork.add_edges_from(
#                     [(best_parent, up_node, network[best_parent][up_node])]
#                 )
#                 up_node = best_parent
#         else:
#             done = True
#     restrictedNetwork.Clean_Up()
#     # Debug: check if we find a valid network
#     restrictedNetwork.IsANetwork()
#     restrictedNetwork.ScoreEdges()
#     return restrictedNetwork
