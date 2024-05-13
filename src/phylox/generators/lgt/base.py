"""
A module for generating (random) LGT phylogenetic networks.

By Joan Carles Pons, Celine Scornavacca, Gabriel Cardona
With their paper: Generation of Level-k LGT Networks
PMID: 30703035 DOI: 10.1109/TCBB.2019.2895344

Adapted by RemieJanssen to output networks with a given number of leaves and reticulations
"""

import networkx as nx
import numpy as np
from networkx.utils.decorators import np_random_state, py_random_state

from phylox import DiNetwork


def _last_node(net):
    """
    The lgt generator names its nodes with integers, so the last node added
    is the node represented with the maximum integer.

    :param net: a phylox DiNetwork
    :return: the node of the network with the highest integer associated to it.
    """
    return max(net.nodes())


def speciate(net, leaf):
    """
    Add a new leaf to the network by a 'speciation event' at leaf `leaf` of the network.
    This is done by adding two outgoing edges to the `leaf` node.

    :param net: a phylox DiNetwork
    :param leaf: leaf of the network `net`
    :return: None
    """
    l = _last_node(net)
    net.add_edge(leaf, l + 1)
    net.add_edge(leaf, l + 2)


def lgt(net, leaf1, leaf2):
    """
    Adds a 'lateral gene transfer' to the network between leaves `leaf1` and 
    `leaf2`. This is done by adding an outgoing edge to each of the leaf 
    nodes and adding an edge between these nodes as well.

    :param net: a phylox DiNetwork
    :param leaf1: a leaf of the network `net`
    :param leaf2: a leaf of the network `net`
    :return: None
    """    
    net.add_edge(leaf1, leaf2)
    l = _last_node(net)
    net.add_edge(leaf1, l + 1)
    net.add_edge(leaf2, l + 2)


def leaves(net):
    """
    Returns a list of leaves of the network.

    :param net: a phylox DiNetwork
    :return: a list of leaf nodes of the network
    """
    return [u for u in net.nodes() if net.out_degree(u) == 0]

def reticulations(net):
    """
    Returns a list of reticulation nodes of the network.

    :param net: a phylox DiNetwork
    :return: a list of reticulation nodes of the network
    """
    return [v for v in net.nodes() if net.in_degree(v) == 2]


def internal_blobs(net):
    """
    Returns a list of internal blobs of the network.
    An internal blob is a biconnected component of 
    the underlying undirected graph of the network, after removing the leaves.

    :param net: a phylox DiNetwork
    :return: a list of internal blobs of the network
    """
    internal_nodes = set([u for u in net.nodes() if net.out_degree(u) > 0])
    blobs = list(nx.biconnected_components(nx.Graph(net)))
    blobs = [bl for bl in blobs if len(bl) > 2]
    nodes_in_blobs = set().union(*blobs)
    nodes_not_in_blobs = internal_nodes - nodes_in_blobs
    blobs.extend([set([u]) for u in nodes_not_in_blobs])
    return blobs


def compute_hash(net):
    """
    Returns a dictionary, mapping each internal node to its blob, 
    and each leaf to the blob directly adjacent to it.

    :param net: a phylox DiNetwork
    :return: a dictionary mapping each leaf to an internal blob.
    """
    mapping_blobs = {}
    blobs = internal_blobs(net)
    for blob in blobs:
        for node in blob:
            mapping_blobs[node] = blob

    mapping = {}
    for l in leaves(net):
        parent = list(net.predecessors(l))[0]
        mapping[l] = mapping_blobs[parent]
    return mapping


def internal_and_external_pairs(net):
    """
    Returns a list internal pairs, and a list of external pairs.
    A pair of leaves is internal if they are adjacent to the same internal blob, 
    the pair is external otherwise.
    
    :param net: a phylox DiNetwork
    :return: a list internal pairs, and a list of external pairs
    """
    lvs = leaves(net)
    pairs = [(l1, l2) for l1 in lvs for l2 in lvs if l1 != l2]
    mapping = compute_hash(net)
    internal_pairs = []
    external_pairs = []
    for pair in pairs:
        if mapping[pair[0]] == mapping[pair[1]]:
            internal_pairs.append(pair)
        else:
            external_pairs.append(pair)
    return internal_pairs, external_pairs


@np_random_state("seed")
def random_pair(net, wint, wext, seed=None):
    """
    Randomly returns a pair of leaves, weighted by two parameters: a pair has
    weight `wint` if it is an internal pair, and `wext` otherwise.
    
    :param net: a phylox DiNetwork
    :param wint: The weight (float) of an internal leaf pair
    :param wext: The weight (float) of an external leaf pair
    :return: a pair of leaves of `net`
    """
    int_pairs, ext_pairs = internal_and_external_pairs(net)
    probabilities = [wint] * len(int_pairs) + [wext] * len(ext_pairs)
    probabilities = np.array(probabilities) / sum(probabilities)
    index = seed.choice(range(len(int_pairs) + len(ext_pairs)), p=probabilities)
    return (int_pairs + ext_pairs)[index]


@np_random_state("seed")
def simulation_1(num_steps, prob_lgt, wint, wext, seed=None):
    """
    Generates a random DiNetwork starting with a tree on two leaves,
    and adding `num_steps` events (either speciation or lgt). The type
    of event is chosen independently at random in each step with 
    probability `prob_lgt` for an lgt event, and probability `1-prob_lgt`
    for a speciation event. 
    
    An lgt is added between a pair of leaves, the pair
    is chosen by weighting a leaf pair with weight `wint` if it is an
    internal pair, and weight `wext` if it is an external pair.
    A pair of leaves is internal if they are adjacent to the same internal blob, 
    the pair is external otherwise.

    :param num_steps: the total number of speciation and lgt events
    :param prob_lgt: the probability of a lgt event (a float in [0,1])
    :param wint: The weight (float) of an internal leaf pair
    :param wext: The weight (float) of an external leaf pair
    :param seed: The random seed
    :return: a phylox DiNetwork
    """

    net = DiNetwork(edges=[(0, 1), (1, 2), (1, 3)])
    net.add_edge(1, 2)
    net.add_edge(1, 3)
    for i in range(num_steps):
        event = seed.choice(["spec", "lgt"], p=[1 - prob_lgt, prob_lgt])
        # event = np.random.choice(['spec','lgt'],p=[1-prob_lgt, prob_lgt])
        if event == "spec":
            l = seed.choice(leaves(net))
            speciate(net, l)
        else:
            pair = random_pair(net, wint, wext, seed=seed)
            lgt(net, pair[0], pair[1])
    return net


@np_random_state("seed")
def generate_network_lgt(leaves_goal, retics_goal, wint, wext, seed=None):
    """
    Generates a random DiNetwork with a given number of leaves and reticulations.
    The number of leaves and reticulations uniquely determines the number of 
    speciation and lgt events, but not their order. Hence, the order of the events
    is determined first (at random), and then the events are all applied randomly.

    An lgt edge is added between a pair of leaves, the pair
    is chosen by weighting a leaf pair with weight `wint` if it is an
    internal pair, and weight `wext` if it is an external pair.
    A pair of leaves is internal if they are adjacent to the same internal blob, 
    the pair is external otherwise.

    :param leaves_goal: number of leaves in the network
    :param retics_goal: number of reticulations in the network
    :param wint: The weight (float) of an internal leaf pair
    :param wext: The weight (float) of an external leaf pair
    :return: a phylox DiNetwork with the given number of leaves and reticulations
    """
    original_leaves_goal = leaves_goal
    # if goal is 1 leaf: pretend we need two leaves, and connect them again later
    if leaves_goal == 1:
        if retics_goal == 0:
            return DiNetwork(edges=[(0, 1)])
        leaves_goal = 2

    # pick a number of extant lineages for each LGT event independently
    retics_at_lineage = dict()
    for r in range(retics_goal):
        lin = seed.choice(range(2, leaves_goal + 1))
        retics_at_lineage[lin] = retics_at_lineage.get(lin, 0) + 1

    # create a network and apply the events
    network = DiNetwork(edges=[(0, 1), (1, 2), (1, 3)])
    for _ in range(retics_at_lineage.get(2, 0)):
        pair = random_pair(network, wint, wext, seed=seed)
        lgt(network, pair[0], pair[1])
    for i in range(3, leaves_goal + 1):
        l = seed.choice(leaves(network))
        speciate(network, l)
        for _ in range(retics_at_lineage.get(i,0)):
            pair = random_pair(network, wint, wext, seed=seed)
            lgt(network, pair[0], pair[1])

    if original_leaves_goal != 1:
        network._set_leaves()
        return network
    # if the goal was 1 leaf, join the two leaves
    unused_node = _last_node(network)
    for leaf in network.leaves:
        leaf_parent = network.parent(leaf)
        network.remove_node(leaf)
        network.add_edge(leaf_parent, unused_node)
    network.add_edge(unused_node, unused_node + 1)
    network._set_leaves()
    return network


@np_random_state("seed")
def generate_network_lgt_conditional(n, k, prob_lgt, wint=1, wext=1, max_tries=1000, seed=None):
    """
    Generates a random DiNetwork with a given number of leaves and reticulations.
    The number of leaves and reticulations uniquely determines the number of 
    speciation and lgt events, but not their order. In contrast to `simulation_3`
    the order of the events is not determined beforhand, but generated like in
    `simulation_1` this guarantees a draw according to simulation_1 conditional on
    the number of leaves and reticulations

    An lgt is added between a pair of leaves, the pair
    is chosen by weighting a leaf pair with weight `wint` if it is an
    internal pair, and weight `wext` if it is an external pair.
    A pair of leaves is internal if they are adjacent to the same internal blob, 
    the pair is external otherwise.

    :param n: number of leaves
    :param k: number of reticulations
    :param alpha: parameter for the weight of internal edges
    :param beta: parameter for the weight of external edges
    :param max_tries: maximum number of tries to generate a network
    :param seed: seed for the random number generator
    :return: a network with the given number of leaves and reticulations
    """

    original_leaves_goal = n
    # if goal is 1 leaf: pretend we need two leaves, and connect them again later
    if n == 1:
        if k == 0:
            return DiNetwork(edges=[(0, 1)])
        n = 2

    for _ in range(max_tries):
        network = simulation_1(n+k-2, prob_lgt, wint, wext, seed=seed)
        if len(reticulations(network)) == k:
            if original_leaves_goal != 1:
                network._set_leaves()
                return network
            unused_node = _last_node(network)
            for leaf in network.leaves:
                leaf_parent = network.parent(leaf)
                network.remove_node(leaf)
                network.add_edge(leaf_parent, unused_node)
            network.add_edge(unused_node, unused_node + 1)
            network._set_leaves()
            return network
            
    raise Exception(
        "Could not generate network with %d leaves and %d reticulations" % (n, k)
    )
