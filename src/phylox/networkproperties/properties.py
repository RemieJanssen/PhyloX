import math

import networkx as nx

from phylox import find_unused_node
from phylox.cherrypicking import is_second_in_reducible_pair


def count_reducible_pairs(network):
    """returns the number of reducible pairs in the network
    split up by number of cherries and number of reticulated cherries.
    """
    cherries = []
    reticulate_cherries = []
    for node in network.nodes:
        if network.out_degree(node) == 0:
            pair = is_second_in_reducible_pair(network, node)
            if pair != False:
                for parent in network.predecessors(pair[0]):
                    if network.out_degree(parent) == 1:
                        reticulate_cherries += [pair]
                        break
                else:
                    cherries += [pair]
    return {
        "cherries": len(cherries) / 2,
        "reticulate_cherries": len(reticulate_cherries),
    }


def blob_properties(network):
    """returns a list of all blobs of the network and their properties.
    Each blob is a pair (blob_size, blob_level) where blob_size is the
    number of nodes in the blob and blob_level is the number of
    reticulations in the blob.
    """
    blob_properties = []
    # For each biconnected component
    for bicomponent in nx.biconnected_components(network.to_undirected()):
        # A bicomponent is a blob if it consists of at least 2 nodes
        if len(bicomponent) > 2:
            #            blob = network.subgraph(bicomponent)
            # count reticulations in blob
            retics = 0
            for node in bicomponent:
                if network.in_degree(node) == 2:
                    retics += 1
            blob_size = len(bicomponent)
            blob_level = retics
            blob_properties += [(blob_size, blob_level)]
    return blob_properties


def level(network):
    """returns the level of the network"""
    blobs = blob_properties(network)
    return max([blob[1] for blob in blobs])


def b2_balance(network, connect_roots=False):
    """returns the B_2 balance of the network"""
    balance = 0

    # Initiate probabilities
    roots = network.roots
    if len(roots) == 0:
        return 0
    if len(roots) > 1:
        if connect_roots:
            new_root = find_unused_node(network)
            for root in roots:
                network.add_edge(new_root, root)
        else:
            raise ValueError("Network has more than one root")
        root = new_root
    else:
        root = list(roots)[0]

    probabilities = dict()
    probabilities[root] = 1
    highest_nodes = []
    for node in network.successors(root):
        if network.in_degree(node) == 1:
            highest_nodes += [node]

    # Calculate the probabilities of reaching each node with a uniform random directed walk
    while highest_nodes:
        current_node = highest_nodes.pop()
        # Calculate the probability of current_node
        prob = 0
        for parent in network.predecessors(current_node):
            prob += (1 / float(network.out_degree(parent))) * probabilities[parent]
        probabilities[current_node] = prob

        # Update the set of highest nodes
        for child in network.successors(current_node):
            for parent in network.predecessors(child):
                if parent not in probabilities:
                    break
            else:
                highest_nodes += [child]

        #
        if network.out_degree(current_node) == 0:
            balance -= prob * math.log(prob, 2)

    return balance
