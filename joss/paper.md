---
title: 'PhyloX: A Python package for complete phylogenetic network workflows'
tags:
  - Python
  - Bio-informatics
  - Phylogenetics
  - Graph theory
authors:
  - name: Remie Janssen
    orcid: 0000-0002-5192-1470
    affiliation: 1
affiliations:
 - name: National Institute for Public Health and the Environment, Bioinformatics and Computing group, Bilthoven, The Netherlands
   index: 1
date: 25 April 2024
bibliography: paper.bib
---

# Summary

PhyloX is a Python package with tools for generating, manipulating, and analyzing phylogenetic networks. It uses the NetworkX package [@SciPyProceedings_11] for basic graph operations. This has the added benefit that the powerful graph tools from NetworkX can be used directly on the phylogenetic networks as well. The aim of the package is to be of general use to phylogenetic network researchers, with a current focus on I/O, random generation of networks, cherry-picking methods, rearrangement operations, and the identification of classes and properties of networks.



# Phylogenetic networks

In the study of the evolutionary history of biological species and languages, it is common to represent putative histories using graphs. Traditionally, at least in biology, these graphs were most often trees, such as the well known tree drawn by Charles Darwin in one of his notebooks. A tree like this called a phylogenetic tree. In some cases, the evolutionary history includes complex processes like horizontal gene transfer and hybridization. These processes cause a *reticulate* (i.e. network-like) structure in the evolutionary history, which requires phylogenetic networks to be used for representing the evolutionary histories.

A *directed phylogenetic network* (e.g., [@huson2010phylogenetic]) is a directed acyclic graph with four types of nodes:
 - a *root*: an in-degree 0, out-degree 1 node,
 - a labelled set of *leaves*: in-degree 1, out-degree 0 nodes,
 - a set of *reticulation nodes*: in-degree $>1$, out-degree 1 nodes,
 - a set of *tree nodes*: in-degree 1, out-degree $>1$ nodes.

A network is *binary* if each reticulation node has in-degree 2, and each tree node has in-degree 2. An *undirected phylogenetic network* is the underlying undirected graph of a directed phylogenetic network, retaining the labelling of the leaf nodes.

## Network properties
When analysing or comparing phylogenetic networks or phylogenetic network methods, it is often useful to extract some (numerical) parameters from the networks. Some of the most used properties are the *reticulation number* (the number of reticulation nodes in a binary network), the number of *blobs* (biconnected components of the network), and the *level* (the maximum reticulation numberamong all blobs of the network). Of course, there are more, like the recently introduced $B_2$-balance index of the network [@franccois2021revisiting].

## Classes of networks
In research on phylogenetic networks, it is common to restrict attention to some well-known classes of phylogenetic networks. These classes put additional restrictions on the definition of a network, for the benefit of computational efficiency, to model certain biological restraints, or for both.

Kong et al. [@kong2022classes] gives a good overview of most well-known classes of directed phylogenetic networks, and their biological interpretation. For example, tree-child networks are networks in where each ancestral species has a descendant among the extant taxa (the leaves) through only mutation in the network [@4407681]. Mathematically, they are characterized as networks in which each non-leaf node has at least one child that is not a reticulation node.

## Cherry-picking
A basic structure in any network or tree is the *cherry*, a pair of leaves with a common parent. A modified version often found in phylogenetic networks is the *reticulated cherry*, an ordered pair of leaves $(x, y)$ which are related through the three edges $(p_x, x)$, $(p_y, y)$, and $(p_y, p_x)$.

A common modification to a phylogenetic network is to *pick* or ( or *reduce*) a cherry or reticulated cherry. To pick a cherry $(x,y)$, one removes the leaf $x$ from the network together with its incoming edge, and then suppresses the resulting degree 2 node if the shared parent of $x$ and $y$ had out-degree 2. *Suppressing* a degree 2 node consists in removing the node and its two incident edges, and replacing them all with one new edge. To pick a reticulated cherry $(x,y)$, one removes the edge $(p_y, p_x)$, and supresses all resulting degree 2 nodes. The reverse action of picking a cherry, is called *adding a cherry*.

These modifications are used in computational tools, for example to reconstruct networks from ancestral profiles [@ERDOS201933; @bai2021defining; @cardona2024comparison], to check whether one tree-child network is contained in another [@janssen2021cherry], or to combine multiple trees into one network [@linz2019attaching; @van2022practical]. Their versatile use also lead to the introduction of the class of *orchard networks* [@janssen2021cherry; @ERDOS201933]. This class contains all networks that can be reduced to a single leaf using cherry-picking operations. Networks from this class can be interpreted as trees with horizontal gene transfer arcs [@van2022orchard].

## Rearranging networks
For phylogenetic inference problems, it is common to use heuristics that search through a space of networks. Such a space of networks takes the shape of a graph, whose objects are all networks with a common set of leaf labels (the sampled taxa), and sometimes also a set number of reticulations. The edges of the graph correspond to small changes made to a network: there is an edge between two networks if one can make a modification to one of the networks to arrive at the second network.

The modifications that are allowed are well-defined as types of *rearrangement* opterations. These operations can be *horizontal*, keeping the reticulation number the same, or *vertical*, changing the reticulation number. Most horizontal operations are a variation on or restriction of the *rooted subtree prune and regraft* (*rSPR*) operation [@bordewich2017lost; @gambette2017rearrangement].

The names for vertical moves have not been standardized, but they generally do the same. A vertical move that removes a reticulation removes an incoming edge of a reticulation node, and then suppresses the resulting degree 2 nodes. A vertical move that adds a reticulation does the reverse: it *subdivides* two edges of the network, and adds a new edge between the two new degree 2 nodes (e.g. [@bordewich2017lost]).

As mentioned, rearrangement moves can be used to traverse a space of networks. This is used, for example, to sample posterior distributions in Bayesian analyses [@wen2016bayesian, @zhang2018bayesian] and to find networks that score high on a maximum likelihood criterium [@wen2016reticulate; @yu2014maximum].

## Generating networks
To test phylogenetic network methods, one either needs to source or create a test set of networks. Creating them is often the simpler option, so methods to randomly generate phylogenetic networks are ready at hand. Moreover, these methods are often based on evolutionary models that are defined on a high level, i.e. with explicit events for  with processes such as speciation, extinction, and hybridization.

The paper [@janssen2021comparing] contains a comparison of several 'generators', including several previously existing ones (e.g., [@pons2019generation] and [@zhang2018bayesian]) and a new extention of a tree generator to networks.

## Representing networks
Because phylogenetic networks are graphs, a common representation is as a list of edges. Another commonly used representation is the extended Newick format [@cardona2008extended]. The extended Newick notation has a further extension (Rich Newick format) that adds numerical parameters to the edges of the network, such as the branch length, and the inheritance probability (for incoming edges of a reticulation node) [@riceRichNewick; @wen2018inferring].



# PhyloX Functionality

PhyloX is equipped to handle all the aspects pf phylogenetic networks mentioned in the previous section. It is written primarily for explorative research into algorithmic aspects of phylogenetic networks, although application focused implementations can also be realized with it. An example is the software [@https://doi.org/10.4121/c679cd3c-0815-4021-a727-bcb8b9174b27.v1] for the paper [@bernardini2023constructing], which uses cherry-picking methods in combination with machine learning to efficiently combine a large number of trees into a phylogenetic network. This software shares some of its basic code with the [cherrypicking module](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.cherrypicking.html) and the [generators module](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.generators.html) of PhyloX.

## I/O
PhyloX handles all stages of a phylogenetic workflow involving networks. This starts and ends with input/output of networks. The [DiNetwork class](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.dinetwork.html), which is used to represent phylogenetic networks in PhyloX, inherits from the [DiGraph](https://networkx.org/documentation/stable/reference/classes/digraph.html) class of NetworkX [@SciPyProceedings_11]. Hence, `phylox.DiNetwork` objects can simply be created using the API of `networkx.DiGraph`, and adding labels to the leaves:
```python
from phylox import DiNetwork
from phylox.constants import LABEL_ATTR

network = DiNetwork()
network.add_edges_from(((0,1),(1,2),(1,3)))
network.nodes[2][LABEL_ATTR] = "leaf1"
network.nodes[3][LABEL_ATTR] = "leaf2"
```

The same can be achieved with a modified initialization of DiNetwork:
```python
from phylox import DiNetwork

network = DiNetwork(
    edges=((0,1),(1,2),(1,3)),
    labels=[(2,"leaf1"), (3,"leaf2")]
)
```

Alternatively, the network can be initialized from a Newick string with

```python
from phylox import DiNetwork

network = DiNetwork.from_newick("((leaf1,leaf2));")
```

For output, it is also possible to use functionality from NetworkX. For example, it is possible to output the list of edges or to create a drawing of the network. Of course, output as Newick string is also available with PhyloX (with `network.newick()` for a network called `network` as in the example code blocks above). This outputs all edge information in rich Newick format by default, but can also be forced to output an extended Newick string without edge information.

## Generating networks
Networks can also be generated randomly in PhyloX, which can be utilized to create test sets for new methods. The implemented generators are based on the code from [@janssen2021comparing]. These include generators based on evolutionary models, such as the [LGT generator](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.generators.lgt.html) and the [ZODS generator](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.generators.zods.html) based on [@pons2019generation] and [@zhang2018bayesian], but also a [Metropolis-Hastings sampler] enabling uniform sampling from classes of networks.

The latter makes use of a large part of the functionality of PhyloX, especially when sampling orchard networks: after generating or choosing a starting network, the [`phylox.generators.mcmc.sample_mcmc_networks`](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.generators.mcmc.html) randomly traverses the space of phylogenetic networks using the [rearrangement module](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.rearrangement.html), and rejects proposals if the resulting network is not orchard using the [cherry-picking module](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.cherrypicking.html).
```python
from phylox.generators.randomTC import generate_network_random_tree_child_sequence
from phylox.generators.mcmc import sample_mcmc_networks
from phylox.classes import is_orchard
from phylox.rearrangement.move import MoveType

# Generate an arbitrary orchard network with 10 leaves and 5 reticulations
start_network = generate_network_random_tree_child_sequence(10, 5, seed=4321)
# Generate 100 orchard networks with 10 leaves and 5 reticulations
sampled_networks = sample_mcmc_networks(
    start_network,
    {MoveType.TAIL: 0.5, MoveType.HEAD: 0.5},
    number_of_samples=100,
    burn_in=5,
    restriction_map=is_orchard,
    add_root_if_necessary=True,
    correct_symmetries=False,
    seed=1234,
)
# Write the sampled networks to a file
with open("sampled_networks.nwk", "w") as f:
    for network in sampled_networks:
        f.write(network.newick() + "\n")
```

For this sampler to work correctly, the space of networks that is sampled from needs to be connected. That is, it has to be possible to transform each network into each other network in the space using the selected rearrangement moves. In the example above, this means that the space of orchard network with 10 leaves and 5 reticulations needs to be connected under tail moves and head moves (i.e. rSPR moves).

This is something the user needs to check or prove themselves, as it is not viable to check this computationally. Fortunately, such connectivity results have been studied in detail [@klawitter2020spaces; @thesis_janssen; @van2022orchard; @ERDOS2021205]. For example, the result needed to prove that this example is correct can be found in [@van2022orchard].

## Comparing networks
Based on all the properties above, PhyloX provides a toolkit to compare networks. For example, it can be used to determine whether two networks are [isomorphic](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.isomorphism.base.is_isomorphic.html#phylox.isomorphism.base.is_isomorphic) (i.e., the same); whether they have the same [properties](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.networkproperties.html): level, number of blobs, reticulation number, and number of (reticulated) cherries; whether one is [contained](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.cherrypicking.tree_child_sequences.html#phylox.cherrypicking.tree_child_sequences.tree_child_network_contains) in the other, if both are tree-child; and whether they are similar with respect to a [rearangement distance](https://phylox.readthedocs.io/en/v1.0.3/_autosummary/phylox.rearrangement.exact_distance.html#module-phylox.rearrangement.exact_distance).



# Statement of Need

Currently, there is no Python package that enables a full workflow for analysing properties and methods of phylogenetic networks. Isolated scripts for this purpose do appear on GitHub or as pseudo-code regularly, most often as part of publications studying one method or one property [@zhang2018bayesian; @pons2019generation; @janssen2020combining; @janssen2020linear; @janssen2021rearranging]. Combining such scripts requires substantial work, for example because the phylogenetic networks themselves are represented by different Python classes with their own methods.

This package, PhyloX, aims to bring these scripts together: it standardizes implementations of several basic objects related to phylogenetic networks, such as the networks themselves, the labelling of the nodes, and rearrangement moves. It currently implements a limited but important set of basic functions: I/O for networks (e.g. lists of edges and extended newick format), network generation for test sets, comparing networks resulting from reconstruction methods, and computing several well-used network properties such as the reticulation number, the level, and the number of cherries.

## Related packages
As mentioned above, there are currently no Python packages that enable a complete workflow for phylogenetic networks. However, some Python packages are available that enable part of this workflow or a very similar one. In this section, we compare the functionality of several of these packages to PhyloX, focussing only on usability for phylogenetic networks.

### PhyloNetwork
Like PhyloX, [PhyloNetwork](https://github.com/bielcardona/PhyloNetwork) is a Python package based on NetworkX. It has a richer implementation for phylogenetic trees than PhyloX. For example, it includes more tree-specific rearrangement moves, the calculation of node properties such as the latest common ancestor (LCA), and some presets for drawing networks.

However, it has very few methods for phylogenetic networks and most of those methods are also included in PhyloX. Another advantage of using PhyloX over PhyloNetwork is the inclusion of explicit random seeds. This is an important factor for the reproducibility of research.

Note that code from PhyloNetwork and PhyloX may be easy to combine, as both use NetworkX to implement the phylogenetic network class.

### Biopython - Phylo
This phylogenetics module, [Phylo](https://biopython.org/wiki/Phylo) [@talevich2012bio], of the Biopython package [@10.1093/bioinformatics/btp163] is built for phylogenetic analyses in Python. However, it is set up for phylogenetic trees only. The encoding of trees as sets of [clades](https://github.com/biopython/biopython/blob/17a9a5e41bafd9a85df18d4210e9293707d9e369/Bio/Phylo/PhyloXML.py#L321) does not easily allow extension to networks, which makes it unsuitable to use for these phylogenetic networks methods.

### DendroPy
Like Biopython's phylogenetics package, the [DendroPy](https://github.com/jeetsukumaran/DendroPy) package focuses on phylogenetic trees [@sukumaran2010dendropy]. Unlike Biopython, the implementation of the trees in DendroPy is graph based, making it more suited for analyses of phylogenetic networks. This could still require large changes, as some properties of trees are built into the code on a fairly fundamental level, such as each node having (at most one) [parent node](https://github.com/jeetsukumaran/DendroPy/blob/cc82ab774ed83831b5c5125278d88c3c614c2d8a/src/dendropy/datamodel/treemodel/_node.py#L55C14-L55C26).



# Availability

The code of PhyloX is available as an open source project on [GitHub](https://github.com/RemieJanssen/PhyloX) under the BSD 3-Clause license. The package is also available via [PyPI](https://pypi.org/project/phylox/), so it can be installed via pip (or pip in conda), and updates to the release branch are automatically converted into new versions of the package. The releases are recorded in [Zenodo](https://zenodo.org/records/10122404) so persistent identifiers can be used to cite specific releases of the software. When citing this software, please make sure to also cite the original source of the code, which is mentioned in the [documentation](https://phylox.readthedocs.io/) of each method or class.

\section{Acknowledgements}
Most of the code has been written in the form of separate scripts during the author's PhD project, which was conducted under Leo van Iersel's \orcidlink{0000-0001-7142-4706} Vidi grant: 639.072.6

Anyone willing to contribute is very welcome to do so via pull requests and issues on GitHub!



# References
