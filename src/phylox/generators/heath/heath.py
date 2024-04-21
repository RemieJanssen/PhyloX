"""
This module contains the code to generate phylogenetic networks using the
Heath model as it was called in the paper:
"Comparing the topology of phylogenetic network generators"
by Remie Janssen and Pengy Liu (2021).

This heath model is an extention of the model by Heath et al for generating trees:
"Taxon sampling affects inferences of macro-evolutionary processes from phylogenetic trees" by
Heath TA, Zwickl DJ, Kim J, Hillis DM (2008).

The main difference is that hybridization and HGT events are allowed, the rates of these events are
handled similarly to the speciation and extinction rates in the original model.
"""

import itertools
import math
import random
import sys

import networkx as nx
import numpy as np
import scipy.stats

from networkx.utils.decorators import np_random_state, py_random_state

from phylox import DiNetwork
from phylox.constants import LENGTH_ATTR, PROBABILITY_ATTR, LABEL_ATTR


def gamma_distribution_pdf(value, mean, shape):
    """
    A reparameterization of the Gamma Distribution from scipy stats,
    with mean parameter instead of a shape parameter.

    """
    scale = mean / shape
    return scipy.stats.gamma.pdf(value, shape, 0, scale)


@np_random_state("seed")
def update_rate(parent_rate, prior_mean, prior_shape, update_shape, seed=None):
    """
    Generates the new value for a rate used in the Heath generator based on
     - the current value,
     - the multiplicative update factor `np.random.gamma(<update_shape>, 1 / <update_shape>)`
     - the prior Gamma Distribution for the rate type `GammaDistributionPDF(rate_value, prior_mean, prior_shape)`
    """

    proposed_rate = parent_rate * seed.gamma(update_shape, 1 / update_shape)
    if seed.random() < gamma_distribution_pdf(
        proposed_rate, prior_mean, prior_shape
    ) / gamma_distribution_pdf(parent_rate, prior_mean, prior_shape):
        return proposed_rate
    return parent_rate

@np_random_state("seed")
def update_all_rates(
    parent_rates,
    update_shape,
    speciation_rate_mean,
    speciation_rate_shape,
    ext_used,
    extinction_rate_mean,
    extinction_rate_shape,
    hgt_used,
    hgt_rate_mean,
    hgt_rate_shape,
    seed=None,
):
    """
    Updates all rates for the Heath generator using the `update_rate` function.
    """
    sp_rate = update_rate(
        parent_rates[0], speciation_rate_mean, speciation_rate_shape, update_shape, seed=seed
    )
    ext_rate = 0
    if ext_used:
        ext_rate = update_rate(
            parent_rates[1], extinction_rate_mean, extinction_rate_shape, update_shape, seed=seed
        )
    hgt_rate = 0
    if hgt_used:
        hgt_rate = update_rate(
            parent_rates[2], hgt_rate_mean, hgt_rate_shape, update_shape, seed=seed
        )
    return (sp_rate, ext_rate, hgt_rate)


def graph_distance_to_hybridization_rate(
    distance,
    hybridization_left_bound,
    hybridization_right_bound,
    hybridization_left_rate,
    hybridization_right_rate,
):
    """
    Returns a hybdridization rate for a given distance.
    The distance represents the evolutionary distance between two extant species.
    The function is used in the Heath generator to determine the hybdridization rate between two current taxa,
    where the distance is a weighthed sum of all up-down distances between the two taxa in the current network.
    The dependency on the distance takes the following shape:
     - hybridization_left_bound:  l
     - hybridization_right_bound: r
     - hybridization_left_rate:   rl
     - hybridization_right_rate:  rr.

    ::

           |
        lr +---
           |   \\
           |    \\
           |     \\
        rr +      -----
           |
         0 +---+----+-----
           0   l    r

    where the distance is on the x-axis, and the hybridization rate on the y-axis.
    """
    if distance <= hybridization_left_bound:
        return hybridization_left_rate
    elif distance >= hybridization_right_bound:
        return hybridization_right_rate
    return hybridization_left_rate + (
        hybridization_right_rate - hybridization_left_rate
    ) * (distance - hybridization_left_bound) / (
        hybridization_right_bound - hybridization_left_bound
    )


def restrict_network_to_leaf_set(network, leaves_to_keep):
    """
    Removes all leaves in the network that are not in the `leaves_to_keep` container.
    Then cleans up the network, by iteratively removing out-degree 0 nodes that are not in the `leaves_to_keep` set, and suppressing in-degree 1 out-degree 1 nodes.
    Modifies the network in place and returns it.
    """
    # find leaves to remove
    leaves_to_keep = set(leaves_to_keep)
    remove_nodes = set()
    for v in network.nodes():
        if network.out_degree(v) == 0 and v not in leaves_to_keep:
            remove_nodes.add(v)

    # remove the sinks that we don't want to keep
    while remove_nodes:
        removed_node = remove_nodes.pop()
        parents = list(network.predecessors(removed_node))
        network.remove_node(removed_node)
        for p in parents:
            if network.out_degree(p) == 0:
                remove_nodes.add(p)

    # suppress degree-2 nodes
    network = suppress_degree_two_nodes(network)
    return network


def suppress_degree_two_nodes(network):
    """
    Suppresses degree 2 (in-degree 1 out-degree 1) nodes of the network in place and
    returns the network. The length of the new edge is the sum of the lengths of the old two edges.
    If the bottom edge had a probability, then this probability is given to the new edge.
    """
    to_remove = []
    to_check = set(list(network.nodes())[:])
    while to_check:
        v = to_check.pop()
        if network.in_degree(v) == 1 and network.out_degree(v) == 1:
            to_remove += [v]
            parent = list(network.predecessors(v))[0]
            child = list(network.successors(v))[0]
            if network.has_edge(parent, child):
                to_check.add(parent)
                to_check.add(child)
            in_edge_length = network[parent][v][LENGTH_ATTR]
            out_edge_length = network[v][child][LENGTH_ATTR]
            out_edge_prob = network[v][child].get(PROBABILITY_ATTR)
            network.remove_edges_from([(parent, v), (v, child)])
            network.add_weighted_edges_from(
                [(parent, child, in_edge_length + out_edge_length)], weight=LENGTH_ATTR
            )
            if out_edge_prob != None:
                network[parent][child][PROBABILITY_ATTR] = out_edge_prob
    network.remove_nodes_from(to_remove)
    return network


# For each choice of reticulation arcs, calculate the distance between all pairs of taxa
# Add up all these distances for each pair, weighed by the probability of this embedded tree
# The hyb rates are e^(-hybridization_rate*sum_for_pair-offset) for each pair.
#####
# Updating distances between pairs can be done smartly:
#  for a speciation event, the distances simply increase by two times the time to speciation, and the new species copies the distances from its sister species
#  for an extinction event, the distances simply increase by two times the time to speciation, and the distances to the extinct species are removed.
#  for a HGT event, only the rates for the receiving taxon have to be updated. They can be computed from the distances to the receiving and the donating easily
#  for a hyb event, all paths go via exactly one of the parent species, so we can use those distances again.
# Updating hybridization rates is done within the generate_heath_network function
@np_random_state("seed")
def generate_heath_network(
    time_limit=1.0,
    taxa_limit=None,
    update_shape=2.0,
    speciation_rate_mean=2.0,
    speciation_rate_shape=2.0,
    ext_used=True,
    count_extinct=False,
    extinction_rate_mean=1.0,
    extinction_rate_shape=1.0,
    hgt_used=False,
    hgt_rate_mean=None,
    hgt_rate_shape=None,
    hgt_inheritance=0.05,
    hyb_used=False,
    hybridization_left_bound=None,
    hybridization_right_bound=None,
    hybridization_left_rate=None,
    hybridization_right_rate=None,
    simple_output=False,
    seed=None,
):
    """
    Runs a speciation-extinction-HGT-hybridization model for the given time (`time_limit`)
    or until a certain number of extant taxa (`taxa_limit`) is reached.
    If all lineages go extinct before the given time is reached, another attempt is made.
    Each extant taxon has its own speciation, HGT, and extinction rates (`rate=1/mean_time_until_next_event`).
    Hybridization rates are `evolutionary distance` dependent, with a function determined by global parameters

    There are prior speciation/extinction/HGT rate distributions:
    gamma distributions with a given mean and a shape parameter for speciation
    (`speciation_rate_mean` and `speciation_rate_shape`)
    and extinction (`extinction_rate_mean` and `extinction_rate_shape`),
    HGT (Horizontal gene transfer) is turned off by default, bu can be turned on by setting `hgt_used=True`.
    This also requires you to set the paramaters for the HGT rate distribution (`hgt_rate_mean` and `hgt_rate_shape`).
    Extinction can be turned off by setting `ext_used=False`.

    If an HGT event happens for a given taxon, another taxon (including itself) is chosen
    uniformly at random to donate genetic material (uniformly distributed contribution in `[0,max_hgt]`
    where `max_hgt` is determined by the `hgt_inheritance`).

    If hybridization is turned on (`hyb_used=True`), the hybdridization rate between two taxa
    is calculated as a function of the distance between those taxa. This distance is a weighthed
    sum of all up-down distances between the two taxa in the current network.
    The dependency on the distance takes the following shape:
     - `hybridization_left_bound`:  l
     - `hybridization_right_bound`: r
     - `hybridization_left_rate`:   rl
     - `hybridization_right_rate`:  rr.

    ::

           |
        lr +---
           |   \\
           |    \\
           |     \\
        rr +      -----
           |
         0 +---+----+-----
           0   l    r

    where the distance is on the x-axis, and the hybridization rate on the y-axis.

    After speciation or hybridization, each rate of the new lineages is set by multiplying the
    (weighted mean) rate of the parent lineage(s) by a gamma-distributed factor with mean 1 and
    a shape parameter (`update_shape`), and then accepting this rate with a probability
    proportional to the prior distribution for this rate.
    This gives an ultrametric network on the extant species.

    The random seed can be set with the `seed` parameter.

    Returns a network without leaf labels, the set of hybrid nodes,
    the set of extant taxa, and the number of extinct taxa.
    """
    # Initiate the network
    nw = DiNetwork()
    nw.add_node(0)
    extant_taxa = set([0])
    current_node = 1
    no_of_extinct = 0
    hybrid_nodes = dict()
    no_of_hybrids = 0

    # Draw initial rates and distances
    current_speciation_rate = seed.gamma(
        speciation_rate_shape, speciation_rate_mean / speciation_rate_shape
    )
    current_extinction_rate = 0.0
    if ext_used:
        current_extinction_rate = seed.gamma(
            extinction_rate_shape, extinction_rate_mean / extinction_rate_shape
        )
    current_hgt_rate = 0.0
    if hgt_used:
        current_hgt_rate = seed.gamma(
            hgt_rate_shape, hgt_rate_mean / hgt_rate_shape
        )
    current_hybridization_rate = float(0)

    # Set the initial leaf rates per leaf
    leaf_rates = dict()
    leaf_rates[0] = (current_speciation_rate, current_extinction_rate, current_hgt_rate)


    # Force the first event to be a speciation
    total_rate = current_speciation_rate
    distances = dict()
    # Pick a time for the first event
    extra_time = seed.exponential(1 / float(total_rate))
    current_time = extra_time

    while len(extant_taxa) > 0 and (
        (not taxa_limit and current_time < time_limit)
        or (taxa_limit and len(extant_taxa) + count_extinct * no_of_extinct < taxa_limit)
    ):
        random_number = seed.random()
        splitting_leaf = None
        extinction_leaf = None
        hgt_donor_leaf = None
        parent_acceptor = None
        hgt_acceptor_leaf = None
        hyb_pair = None
        if random_number < current_speciation_rate / total_rate:
            ######################
            #     Speciation     #
            ######################
            if not simple_output:
                print("speciation")
            random_number = seed.random() * current_speciation_rate
            for leaf, rates in leaf_rates.items():
                if random_number < rates[0]:
                    splitting_leaf = leaf
                    break
                random_number -= rates[0]
            if splitting_leaf == None:
                if not simple_output:
                    print("error, speciation rate computed wrong")
            nw.add_weighted_edges_from(
                [
                    (splitting_leaf, current_node, 0),
                    (splitting_leaf, current_node + 1, 0),
                ],
                weight=LENGTH_ATTR,
            )
            # Update the rates and distances
            # rates
            leaf_rates[current_node] = update_all_rates(
                leaf_rates[splitting_leaf],
                update_shape,
                speciation_rate_mean,
                speciation_rate_shape,
                ext_used,
                extinction_rate_mean,
                extinction_rate_shape,
                hgt_used,
                hgt_rate_mean,
                hgt_rate_shape,
                seed=seed,
            )
            leaf_rates[current_node + 1] = update_all_rates(
                leaf_rates[splitting_leaf],
                update_shape,
                speciation_rate_mean,
                speciation_rate_shape,
                ext_used,
                extinction_rate_mean,
                extinction_rate_shape,
                hgt_used,
                hgt_rate_mean,
                hgt_rate_shape,
                seed=seed,
            )
            current_speciation_rate += (
                leaf_rates[current_node][0]
                + leaf_rates[current_node + 1][0]
                - leaf_rates[splitting_leaf][0]
            )
            current_extinction_rate += (
                leaf_rates[current_node][1]
                + leaf_rates[current_node + 1][1]
                - leaf_rates[splitting_leaf][1]
            )
            current_hgt_rate += (
                leaf_rates[current_node][2]
                + leaf_rates[current_node + 1][2]
                - leaf_rates[splitting_leaf][2]
            )
            # distances
            if hyb_used:
                for l in extant_taxa:
                    if l != splitting_leaf:
                        pair = (splitting_leaf, l)
                        if pair in distances.keys():
                            new_distance = distances[pair]
                        else:
                            pair = (l, splitting_leaf)
                            new_distance = distances[pair]
                        distances[(l, current_node)] = new_distance
                        distances[(l, current_node + 1)] = new_distance
                        del distances[pair]
                distances[(current_node, current_node + 1)] = 0

            extant_taxa.add(current_node)
            extant_taxa.add(current_node + 1)
            extant_taxa.remove(splitting_leaf)
            del leaf_rates[splitting_leaf]
            current_node += 2

        elif (
            random_number
            < (current_extinction_rate + current_speciation_rate) / total_rate
        ):
            ######################
            #     Extinction     #
            ######################
            if not simple_output:
                print("extinction")

            random_number = seed.random() * current_extinction_rate
            for leaf, rates in leaf_rates.items():
                if random_number < rates[1]:
                    extinction_leaf = leaf
                    break
                random_number -= rates[1]
            if extinction_leaf == None:
                if not simple_output:
                    print("ouch, extinction rate computed wrong")

            # Update the rates and distances
            # rates
            current_speciation_rate -= leaf_rates[extinction_leaf][0]
            current_extinction_rate -= leaf_rates[extinction_leaf][1]
            current_hgt_rate -= leaf_rates[extinction_leaf][2]
            # distances
            if hyb_used:
                for l in extant_taxa:
                    if l != extinction_leaf:
                        if (extinction_leaf, l) in distances.keys():
                            del distances[(extinction_leaf, l)]
                        else:
                            del distances[(l, extinction_leaf)]

            del leaf_rates[extinction_leaf]
            extant_taxa.remove(extinction_leaf)
            no_of_extinct += 1

        elif (
            random_number
            < (current_extinction_rate + current_speciation_rate + current_hgt_rate)
            / total_rate
        ):
            ######################
            #      HGT event     #
            ######################
            if not simple_output:
                print("HGT")

            random_number = seed.random() * current_hgt_rate
            for leaf, rates in leaf_rates.items():
                if random_number < rates[2]:
                    hgt_acceptor_leaf = leaf
                    break
                random_number -= rates[2]
            if hgt_acceptor_leaf == None:
                if not simple_output:
                    print("ouch, hgt rate computed wrong")
            if len(extant_taxa) > 1:
                hgt_donor_leaf = seed.choice(list(extant_taxa - set([hgt_acceptor_leaf])))
                for p in nw.predecessors(hgt_acceptor_leaf):
                    parent_acceptor = p
                nw.add_weighted_edges_from(
                    [
                        (hgt_donor_leaf, hgt_acceptor_leaf, 0),
                        (hgt_donor_leaf, current_node, 0),
                        (hgt_acceptor_leaf, current_node + 1, 0),
                    ],
                    weight=LENGTH_ATTR,
                )
                prob = hgt_inheritance * seed.random()
                nw[parent_acceptor][hgt_acceptor_leaf][PROBABILITY_ATTR] = 1 - prob
                nw[hgt_donor_leaf][hgt_acceptor_leaf][PROBABILITY_ATTR] = prob

                hybrid_nodes[hgt_acceptor_leaf] = no_of_hybrids
                no_of_hybrids += 1
                # Update the rates and distances
                # rates
                leaf_rates[current_node + 1] = update_all_rates(
                    tuple(
                        prob * x + (1 - prob) * y
                        for x, y in zip(
                            leaf_rates[hgt_acceptor_leaf], leaf_rates[hgt_donor_leaf]
                        )
                    ),
                    update_shape,
                    speciation_rate_mean,
                    speciation_rate_shape,
                    ext_used,
                    extinction_rate_mean,
                    extinction_rate_shape,
                    hgt_used,
                    hgt_rate_mean,
                    hgt_rate_shape,
                    seed=seed,
                )
                leaf_rates[current_node] = leaf_rates[hgt_donor_leaf]
                current_speciation_rate += (
                    leaf_rates[current_node + 1][0] - leaf_rates[hgt_acceptor_leaf][0]
                )
                current_extinction_rate += (
                    leaf_rates[current_node + 1][1] - leaf_rates[hgt_acceptor_leaf][1]
                )
                current_hgt_rate += (
                    leaf_rates[current_node + 1][2] - leaf_rates[hgt_acceptor_leaf][2]
                )
                del leaf_rates[hgt_donor_leaf]
                del leaf_rates[hgt_acceptor_leaf]
                # distances
                if hyb_used:
                    for l in extant_taxa:
                        if l != hgt_acceptor_leaf:
                            acceptor_pair = (hgt_acceptor_leaf, l)
                            if acceptor_pair in distances.keys():
                                acceptor_distance = distances[acceptor_pair]
                            else:
                                acceptor_pair = (l, hgt_acceptor_leaf)
                                acceptor_distance = distances[acceptor_pair]
                            donor_pair = (l, hgt_donor_leaf)
                            donor_distance = 0.0
                            if l != hgt_donor_leaf:
                                if donor_pair in distances.keys():
                                    donor_distance = distances[donor_pair]
                                else:
                                    donor_pair = (hgt_donor_leaf, l)
                                    donor_distance = distances[donor_pair]
                                distances[(l, current_node)] = donor_distance
                                del distances[donor_pair]
                            distances[(l, current_node + 1)] = (
                                1 - prob
                            ) * acceptor_distance + prob * donor_distance
                            del distances[acceptor_pair]
                    distances[(current_node, current_node + 1)] = distances[
                        (hgt_donor_leaf, current_node + 1)
                    ]
                    del distances[(hgt_donor_leaf, current_node + 1)]
                extant_taxa.remove(hgt_donor_leaf)
                extant_taxa.remove(hgt_acceptor_leaf)
                extant_taxa.add(current_node)
                extant_taxa.add(current_node + 1)
                current_node += 2

        #        else:
        #            print("trying HGT with only one leaf")
        #             #Do nothing, there is only one leaf.

        else:
            ######################
            #    Hybridization   #
            ######################
            if not simple_output:
                print("hybridization")
            # i.e.: pick two leaf nodes, create a hybrid between these two leaves
            random_number = seed.random() * current_hybridization_rate
            for pair, distance in distances.items():
                pair_rate = graph_distance_to_hybridization_rate(
                    distance,
                    hybridization_left_bound,
                    hybridization_right_bound,
                    hybridization_left_rate,
                    hybridization_right_rate,
                )
                if random_number < pair_rate:
                    hyb_pair = pair
                    break
                random_number -= pair_rate
            if hyb_pair == None and not simple_output:
                if len(extant_taxa) == 1:
                    print("ah, no leaves for hyb")
                else:
                    print("ouch, hybridization rate computed wrong")
            nw.add_weighted_edges_from(
                [
                    (current_node, current_node + 1, 0),
                    (hyb_pair[0], current_node, 0),
                    (hyb_pair[1], current_node, 0),
                    (hyb_pair[0], current_node + 2, 0),
                    (hyb_pair[1], current_node + 3, 0),
                ],
                weight=LENGTH_ATTR,
            )
            prob = seed.random()
            nw[hyb_pair[0]][current_node][PROBABILITY_ATTR] = prob
            nw[hyb_pair[1]][current_node][PROBABILITY_ATTR] = 1 - prob
            hybrid_nodes[current_node] = no_of_hybrids
            no_of_hybrids += 1

            # Update the rates and distances
            # rates
            leaf_rates[current_node + 1] = update_all_rates(
                tuple(
                    prob * x + (1 - prob) * y
                    for x, y in zip(leaf_rates[hyb_pair[0]], leaf_rates[hyb_pair[1]])
                ),
                update_shape,
                speciation_rate_mean,
                speciation_rate_shape,
                ext_used,
                extinction_rate_mean,
                extinction_rate_shape,
                hgt_used,
                hgt_rate_mean,
                hgt_rate_shape,
                seed=seed,
            )
            leaf_rates[current_node + 2] = leaf_rates[hyb_pair[0]]
            leaf_rates[current_node + 3] = leaf_rates[hyb_pair[1]]
            current_speciation_rate += leaf_rates[current_node + 1][0]
            current_extinction_rate += leaf_rates[current_node + 1][1]
            current_hgt_rate += leaf_rates[current_node + 1][2]
            # distances
            # TODO FIXED?!: The order may still be wrong. It seems I already delete some distances when I still need them later.
            # These are probably the distances related to the hybrid parent species.
            for l in extant_taxa:
                if l == hyb_pair[0]:
                    pair_0_distance = 0
                else:
                    pair_0 = (hyb_pair[0], l)
                    if pair_0 not in distances.keys():
                        pair_0 = (l, hyb_pair[0])
                    pair_0_distance = distances[pair_0]
                    distances[(l, current_node + 2)] = pair_0_distance
                #                del distances[pair_0]
                if l == hyb_pair[1]:
                    pair_1_distance = 0
                else:
                    pair_1 = (hyb_pair[1], l)
                    if pair_1 not in distances:
                        pair_1 = (l, hyb_pair[1])
                    pair_1_distance = distances[pair_1]
                    distances[(l, current_node + 3)] = pair_1_distance
                #                del distances[pair_1]
                distances[(l, current_node + 1)] = (
                    prob * pair_0_distance + (1 - prob) * pair_1_distance
                )

            if (hyb_pair[0], hyb_pair[1]) in distances:
                distances[(current_node + 2, current_node + 3)] = distances[
                    (hyb_pair[0], hyb_pair[1])
                ]
            else:
                distances[(current_node + 2, current_node + 3)] = distances[
                    (hyb_pair[1], hyb_pair[0])
                ]
            distances[(current_node + 1, current_node + 2)] = distances[
                (hyb_pair[0], current_node + 1)
            ]
            distances[(current_node + 1, current_node + 3)] = distances[
                (hyb_pair[1], current_node + 1)
            ]

            remove_pairs = []
            for pair in distances:
                if hyb_pair[0] in pair or hyb_pair[1] in pair:
                    remove_pairs += [pair]
            for pair in remove_pairs:
                del distances[pair]

            del leaf_rates[hyb_pair[0]]
            del leaf_rates[hyb_pair[1]]
            extant_taxa.remove(hyb_pair[0])
            extant_taxa.remove(hyb_pair[1])
            extant_taxa.add(current_node + 1)
            extant_taxa.add(current_node + 2)
            extant_taxa.add(current_node + 3)
            current_node += 4

        # Now extend all pendant edges of extant taxa
        if len(extant_taxa) == 0:
            break
        for l in extant_taxa:
            pl = -1
            for p in nw.predecessors(l):
                pl = p
            nw[pl][l][LENGTH_ATTR] += extra_time

        # Compute the new rates
        current_hybridization_rate = 0
        if hyb_used:
            for pair, distance in distances.items():
                distances[pair] += 2 * extra_time
                current_hybridization_rate += graph_distance_to_hybridization_rate(
                    distances[pair],
                    hybridization_left_bound,
                    hybridization_right_bound,
                    hybridization_left_rate,
                    hybridization_right_rate,
                )

        total_rate = (
            current_speciation_rate
            + current_extinction_rate
            + current_hgt_rate
            + current_hybridization_rate
        )

        # Compute the time of the next event
        extra_time = seed.exponential(1 / total_rate)
        current_time += extra_time

    # The following corrects for overshooting the time limit
    extra_time += time_limit - current_time
    # nothing has happened yet, and there is only one node
    if len(nw) == 1:
        nw.add_weighted_edges_from([(0, 1, time_limit)], weight=LENGTH_ATTR)
        extant_taxa = set([1])
    # each leaf has a parent node, and we can extend each parent edge to time_limit
    else:
        for l in extant_taxa:
            pl = -1
            for p in nw.predecessors(l):
                pl = p
            nw[pl][l][LENGTH_ATTR] += extra_time

    return nw, hybrid_nodes, extant_taxa, no_of_extinct
