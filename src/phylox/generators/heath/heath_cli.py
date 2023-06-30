import argparse
import itertools
import math
import random
import sys

import networkx as nx
import numpy as np
import scipy.stats

# from AddEdgesToTree import *
from phylox.generators.heath.heath import RestrictToLeafSet, generate_heath_network


def parse_args():
    parser = argparse.ArgumentParser(
        description="""
        Runs a speciation-extinction-HGT-hybridization model for the given time (--time) or until a certain number of extant taxa (--taxa) is reached. If all lineages go extinct before the given time is reached, another attempt is made.
        
        Each extant taxon has its own speciation, HGT, and extinction rates (1/mean_time_until_next_event). 
        
        The hybridization rate of a pair of species is a function of the weighed distance between these species (sum of all up-down distances, weighed by their probability). 
        
        There are prior speciation/extinction/HGT rate distributions: gamma distributions with a given mean and a shape parameter for speciation (--speciation_parameters) and extinction (--extinction_parameters), HGT is turned off by default. 
        
        If a HGT event happens for a given taxon, another taxon (including itself) is chosen uniformly at random to donate genetic material (uniformly distributed contribution in [0,max_hgt] where max_hgt is determined by the --hgt_inheritance parameter). After speciation or hybridization, each rate of the new lineages is set by multiplying the (weighted mean) rate of the parent lineage(s) by a gamma-distributed factor with mean 1 and a shape parameter (--update-shape-parameter), and then accepting this rate with a probability proportional to the prior distribution for this rate. This gives an ultrametric network on the extant species.
        
        Optional arguments:
            -ti or --time followed by the total length (float) for the network. 
            -ta or --taxa followed by the number of taxa at which the simulation stops. If all lineages go extinct before the given number of taxa is reached, another attempt is made.
            -ce or --count_extinct to also count the extinct taxa as part of the taxa limit.-oe or --only_extant to return the network restricted to the extant leaves
        Rates:
            -sp or --speciation_parameters followed by a mean (float) and a shape parameter (float) for the gamma distribution of the speciation rate.
            -ext or --extinction_parameters followed by a mean (float) and a shape parameter (float) for the gamma distribution of the extinction rate.
            -noext or --no_extinction to turn off extinction altogether.
            -hgt or --hgt_parameters followed by a mean (float) and a shape parameter (float) for the gamma distribution of the HGT rate.
            -upd or --update_shape_parameter followed by a shape parameter (float) for the update gamma distribution.
            -hyb or --hyb_factor followed by four floats for the piecewise linear dependence of hybridization rate on the distance:

        left bound:  l
        right bound: r
        left rate:   rl
        right rate:  rr.
        
        |
        lr+---
        |     \
        |      \
        |       \
        rr+      -----
        |
        0 +---+----+-----
        0   l    r""",
        formatter_class=argparse.RawTextHelpFormatter,
    )

    parser.add_argument(
        "-ti",
        "--time",
        dest="time_limit",
        type=float,
        help="Total length for the network",
    )
    parser.add_argument(
        "-ta",
        "--taxa",
        dest="taxa_limit",
        type=float,
        help="Number of taxa at which the simulation stops",
    )
    parser.add_argument(
        "-ce",
        "--count_extinct",
        action="store_true",
        help="Count extinct taxa as part of the taxa limit",
    )
    parser.add_argument(
        "-sp",
        "--speciation_parameters",
        nargs=2,
        type=float,
        metavar=("mean", "shape"),
        help="Mean and shape parameter for the gamma distribution of the speciation rate",
    )
    parser.add_argument(
        "-ext",
        "--extinction_parameters",
        nargs=2,
        type=float,
        metavar=("mean", "shape"),
        help="Mean and shape parameter for the gamma distribution of the extinction rate",
    )
    parser.add_argument(
        "-noext",
        "--no_extinction",
        action="store_true",
        help="Turn off extinction altogether",
    )
    parser.add_argument(
        "-hgt",
        "--hgt_parameters",
        nargs=2,
        type=float,
        metavar=("mean", "shape"),
        help="Mean and shape parameter for the gamma distribution of the HGT rate",
    )
    parser.add_argument(
        "--hgt_inheritance",
        type=float,
        help="Maximum contribution of a donor taxon to the genome of the recipient taxon",
    )
    parser.add_argument(
        "-hyb",
        "--hybridization_factor",
        nargs=4,
        type=float,
        metavar=("l", "r", "rl", "rr"),
        help="Piecewise linear dependence of hybridization rate on the distance",
    )
    parser.add_argument(
        "-upd",
        "--update_shape_parameter",
        type=float,
        help="Shape parameter for the update gamma distribution",
    )
    parser.add_argument(
        "-s", "--simple", action="store_true", help="Enable simple output"
    )

    args = parser.parse_args()

    params = {
        "time_limit": args.time_limit,
        "taxa_limit": args.taxa_limit,
        "count_extinct": args.count_extinct,
        "speciation_rate_mean": args.speciation_parameters[0]
        if args.speciation_parameters
        else None,
        "speciation_rate_shape": args.speciation_parameters[1]
        if args.speciation_parameters
        else None,
        "extinction_rate_mean": args.extinction_parameters[0]
        if args.extinction_parameters
        else None,
        "extinction_rate_shape": args.extinction_parameters[1]
        if args.extinction_parameters
        else None,
        "ext_used": not args.no_extinction,
        "hgt_rate_mean": args.hgt_parameters[0] if args.hgt_parameters else None,
        "hgt_rate_shape": args.hgt_parameters[1] if args.hgt_parameters else None,
        "hgt_used": args.hgt_parameters is not None,
        "hybridization_left_bound": args.hybridization_factor[0]
        if args.hybridization_factor
        else None,
        "hgt_inheritance": args.hgt_inheritance,
        "hybridization_right_bound": args.hybridization_factor[1]
        if args.hybridization_factor
        else None,
        "hybridization_left_rate": args.hybridization_factor[2]
        if args.hybridization_factor
        else None,
        "hybridization_right_rate": args.hybridization_factor[3]
        if args.hybridization_factor
        else None,
        "hyb_used": args.hybridization_factor is not None,
        "update_shape": args.update_shape_parameter,
        "simple_output": args.simple,
    }
    return params


def main():
    params = parse_args()
    # Find a network
    if params["taxa_limit"]:
        leaves = []
        no_of_extinct = 0
        while (
            len(leaves) + params["count_extinct"] * no_of_extinct
            != params["taxa_limit"]
        ):
            if not params["simple_output"]:
                print("starting over")
            network, hybrid_nodes, leaves, no_of_extinct = generate_heath_network(
                **params
            )
    else:
        print(params)
        network, hybrid_nodes, leaves, no_of_extinct = generate_heath_network(**params)

    # Restrict to extant leaves if wanted
    if only_extant:
        network = RestrictToLeafSet(network, leaves)

    for e in network.edges:
        info = ""
        if e[0] in hybrid_nodes:
            info += "H" + str(hybrid_nodes[e[0]])
        else:
            info += str(e[0])
        info += "  "
        if e[1] in hybrid_nodes:
            info += "H" + str(hybrid_nodes[e[1]])
        else:
            info += str(e[1])
        print(info)  # ,network[e[0]][e[1]]['length'])


if __name__ == "__main__":
    main()
