"""
A module for generating (random) LGT phylogenetic networks.

By Joan Carles Pons, Celine Scornavacca, Gabriel Cardona
With their paper: Generation of Level- k LGT Networks
PMID: 30703035 DOI: 10.1109/TCBB.2019.2895344

Adapted by RemieJanssen to output networks with a given number of leaves and reticulations
"""

import random

import networkx as nx
import numpy as np

from phylox import DiNetwork


def last_node(net):
    return max(net.nodes())


def speciate(net, leaf):
    l = last_node(net)
    net.add_edge(leaf, l + 1)
    net.add_edge(leaf, l + 2)


def lgt(net, leaf1, leaf2):
    net.add_edge(leaf1, leaf2, secondary=True)
    l = last_node(net)
    net.add_edge(leaf1, l + 1)
    net.add_edge(leaf2, l + 2)


def leaves(net):
    return [u for u in net.nodes() if net.out_degree(u) == 0]


def internal_blobs(net):
    internal_nodes = set([u for u in net.nodes() if net.out_degree(u) > 0])
    blobs = list(nx.biconnected_components(nx.Graph(net)))
    blobs = [bl for bl in blobs if len(bl) > 2]
    nodes_in_blobs = set().union(*blobs)
    nodes_not_in_blobs = internal_nodes - nodes_in_blobs
    blobs.extend([set([u]) for u in nodes_not_in_blobs])
    return blobs


def compute_hash(net):
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


def random_pair(net, wint, wext):
    int_pairs, ext_pairs = internal_and_external_pairs(net)
    return random.choices(
        int_pairs + ext_pairs, weights=[wint] * len(int_pairs) + [wext] * len(ext_pairs)
    )[0]


def simulation_1(num_steps, prob_lgt, wint, wext):
    net = nx.DiGraph()
    net.add_edge(1, 2)
    net.add_edge(1, 3)
    for i in range(num_steps):
        event = random.choices(["spec", "lgt"], [1 - prob_lgt, prob_lgt])[0]
        # event = np.random.choice(['spec','lgt'],p=[1-prob_lgt, prob_lgt])
        if event == "spec":
            l = random.choice(leaves(net))
            speciate(net, l)
        else:
            pair = random_pair(net, wint, wext)
            lgt(net, pair[0], pair[1])
    return net


def simulation_3(leaves_goal, retics_goal, wint, wext):
    """
    Simulation 3 for generating networks with a given number of leaves and reticulations
    :param leaves_goal: number of leaves in the network
    :param retics_goal: number of reticulations in the network
    :param wint: weight of internal edges
    :param wext: weight of external edges
    :return: a network with the given number of leaves and reticulations
    """
    original_leaves_goal = leaves_goal
    if leaves_goal == 1:
        if retics_goal == 0:
            return DiNetwork(edges=[(0, 1)])
        # pretend we need two leaves, and connect them again later
        leaves_goal = 2
    # pick a number of extant lineages for each LGT event independently
    retics_at_lineage = dict()
    for r in range(retics_goal):
        lin = random.choice(range(2, leaves_goal + 1))
        if lin in retics_at_lineage:
            retics_at_lineage[lin] += 1
        else:
            retics_at_lineage[lin] = 1
    network = DiNetwork(edges=[(0, 1), (1, 2), (1, 3)])
    if 2 in retics_at_lineage:
        for j in range(retics_at_lineage[2]):
            pair = random_pair(network, wint, wext)
            lgt(network, pair[0], pair[1])
    for i in range(3, leaves_goal + 1):
        l = random.choice(leaves(network))
        speciate(network, l)
        if i in retics_at_lineage:
            for j in range(retics_at_lineage[i]):
                pair = random_pair(network, wint, wext)
                lgt(network, pair[0], pair[1])

    if original_leaves_goal == 1:
        # connect the two leaves
        unused_node = last_node(network)
        new_leaf_edge = (unused_node, unused_node + 1)
        for leaf in network.leaves:
            leaf_parent = network.parent(leaf)
            network.remove_node(leaf)
            network.add_edge(leaf_parent, unused_node)
        network._set_leaves()
    return network


def reticulations(G):
    return [v for v in G.nodes() if G.in_degree(v) == 2]


def generate_network_lgt(n, k, wint=1, wext=1, max_tries=1000):
    """
    Generate a network with a given number of leaves and reticulations
    :param n: number of leaves
    :param k: number of reticulations
    :param alpha: parameter for the weight of internal edges
    :param beta: parameter for the weight of external edges
    :param max_tries: maximum number of tries to generate a network
    :return: a network with the given number of leaves and reticulations
    """

    for i in range(max_tries):
        # naive method for exact number of reticulations,
        # (original by Pons et al., not made for exact number of reticulations)
        # resG = simulation_1(n+k, alpha, 1, beta)
        # faster method:
        network = simulation_3(n, k, wint, wext)
        if len(reticulations(network)) == k:
            return network
    raise Exception(
        "Could not generate network with %d leaves and %d reticulations" % (n, k)
    )
