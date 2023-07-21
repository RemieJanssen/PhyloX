"""
A module for generating (random) phylogenetic networks.

This implements the algorithm from the paper:
Zhang, C., Ogilvie, H.A., Drummond, A.J., Stadler, T.: Bayesian inference of species networks from multilocus
sequence data. Molecular biology and evolution 35(2), 504–517 (2018)

"""

import random

import numpy as np

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR, LENGTH_ATTR


def multree_to_dinetwork(multree, hybrid_nodes):
    """
    Converts a multree to a DiNetwork.

    Parameters
    ----------
    multree : networkx.DiGraph
        The multree to convert.
        Edges must have a LENGTH_ATTR attribute.
    hybrid_nodes : dict
        A dictionary mapping hybrid nodes to their hybrid number.

    Returns
    -------
    network : DiNetwork
        The converted network.
    """

    edges = []
    for e in multree.edges(data=True):
        node_0 = f"#H{hybrid_nodes[e[0]]}" if e[0] in hybrid_nodes else e[0]
        node_1 = f"#H{hybrid_nodes[e[1]]}" if e[1] in hybrid_nodes else e[1]
        edges.append((node_0, node_1, {LENGTH_ATTR: e[2][LENGTH_ATTR]}))
    network = DiNetwork(
        edges=edges,
    )
    return network


def generate_network_zods(time_limit, speciation_rate, hybridization_rate):
    """
    Generates a network with the ZODS model.

    Parameters
    ----------
    time_limit : float
        The time limit of the network.
    speciation_rate : float
        The speciation rate of each lineage in the network.
    hybridization_rate : float
        The hybridization rate of each pair of lineages in the network.

    Returns
    -------
    network : DiNetwork
        The generated network.
    """

    multree = DiNetwork()
    leaves = set([0])
    current_node = 1

    extra_time = np.random.exponential(1 / float(speciation_rate))
    current_time = extra_time
    current_speciation_rate = float(speciation_rate)
    current_hybridization_rate = float(0)
    rate = current_speciation_rate + current_hybridization_rate

    # First create a MUL-tree
    hybrid_nodes = dict()
    no_of_hybrids = 0

    while current_time < time_limit:
        if random.random() < current_speciation_rate / rate:
            # Speciate
            splitting_leaf = random.choice(list(leaves))
            multree.add_weighted_edges_from(
                [
                    (splitting_leaf, current_node, 0),
                    (splitting_leaf, current_node + 1, 0),
                ],
                weight=LENGTH_ATTR,
            )
            leaves.remove(splitting_leaf)
            leaves.add(current_node)
            leaves.add(current_node + 1)
            current_node += 2
        else:
            # Hybridize
            #  i.e.: pick two leaf nodes, merge those, and add a new leaf below this hybrid node.
            merging = random.sample(sorted(leaves), 2)
            l0 = merging[0]
            l1 = merging[1]
            pl0 = -1
            for p in multree.predecessors(l0):
                pl0 = p
            pl1 = -1
            for p in multree.predecessors(l1):
                pl1 = p
            # If pl0==pl1, the new hybridization results in parallel edges.
            if pl0 != pl1:
                no_of_hybrids += 1
                multree.add_weighted_edges_from(
                    [(l0, current_node, 0)], weight=LENGTH_ATTR
                )
                leaves.remove(l0)
                leaves.remove(l1)
                leaves.add(current_node)
                hybrid_nodes[l0] = no_of_hybrids
                hybrid_nodes[l1] = no_of_hybrids
                current_node += 1
        # Now extend all pendant edges
        for l in leaves:
            pl = -1
            for p in multree.predecessors(l):
                pl = p
            multree[pl][l][LENGTH_ATTR] += extra_time
        no_of_leaves = len(leaves)
        current_speciation_rate = float(speciation_rate * no_of_leaves)
        current_hybridization_rate = float(
            hybridization_rate * (no_of_leaves * (no_of_leaves - 1)) / 2
        )
        rate = current_speciation_rate + current_hybridization_rate
        extra_time = np.random.exponential(1 / rate)
        current_time += extra_time

    # Correct for the fact that we might have gone over the time limit
    extra_time -= current_time - time_limit
    # nothing has happened yet, and there is only one node
    if len(multree) == 0:
        multree.add_weighted_edges_from([(0, 1, time_limit)], weight=LENGTH_ATTR)
        leaves = set([1])
    # each leaf has a parent node, and we can extend each parent edge to time_limit
    else:
        multree.add_weighted_edges_from([(-1, 0, 0)], weight=LENGTH_ATTR)
        for l in leaves:
            pl = -1
            for p in multree.predecessors(l):
                pl = p
            multree[pl][l][LENGTH_ATTR] += extra_time

    # Now we have a MUL-tree, but we need to make it a network.
    return multree_to_dinetwork(multree, hybrid_nodes)
