import math

import networkx as nx

from phylox.base import find_unused_node
from phylox.cherrypicking import is_second_in_reducible_pair


def count_reducible_pairs(network):
    """
    finds the number of reducible pairs in the network
    split up by number of cherries and number of reticulated cherries.

    :param network: a phylogenetic network.
    :return: a dictionary with the number of reducible pairs in the network
        keys are "cherries" and "reticulate_cherries"

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.networkproperties.properties import count_reducible_pairs
    >>> network = DiNetwork(
    ...     edges=[(-1,0),(0,1),(1,2),(1,3),(2,3),(2,4),(3,5),(0,6),(6,7),(6,8)],
    ... )
    >>> counts = count_reducible_pairs(network)
    >>> counts["cherries"] == 1 and counts["reticulate_cherries"] == 1
    True
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
    """finds a list of all blobs of the network and their properties.
    Each blob is a pair (blob_size, blob_level) where blob_size is the
    number of nodes in the blob and blob_level is the number of
    reticulations in the blob.

    :param network: a phylogenetic network.
    :return: a list of pairs (blob_size, blob_level)

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.networkproperties.properties import blob_properties
    >>> network = DiNetwork(
    ...     edges=[(1,2),(2,3),(2,4),(3,4),(3,5),(4,6),(6,7),(6,8),(7,8),(7,9),(8,10)],
    ... )
    >>> blob_properties(network)
    [(3, 1), (3, 1)]

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,4),(3,5),(2,5),(3,4),(4,6),(5,7)],
    ... )
    >>> blob_properties(network)
    [(5, 2)]
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
    """
    returns the level of the network

    :param network: a phylogenetic network.
    :return: the level of the network

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.networkproperties.properties import level
    >>> network = DiNetwork(
    ...     edges=[(1,2),(2,3),(2,4),(3,4),(3,5),(4,6),(6,7),(6,8),(7,8),(7,9),(8,10)],
    ... )
    >>> level(network)
    1

    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,4),(3,5),(2,5),(3,4),(4,6),(5,7)],
    ... )
    >>> level(network)
    2
    """

    blobs = blob_properties(network)
    return max([blob[1] for blob in blobs])


def b2_balance(network, connect_roots=False):
    """returns the B_2 balance of the network

    :param network: a phylogenetic network.
    :param connect_roots: if True, connects all roots to a new root.
    :return: the B_2 balance of the network

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.networkproperties.properties import b2_balance
    >>> network = DiNetwork(
    ...     edges=[(0, 1), (1, 2), (1, 3), (2, 4), (2, 5), (3, 6), (3, 7)],
    ... )
    >>> b2_balance(network) == 2
    True
    """
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
