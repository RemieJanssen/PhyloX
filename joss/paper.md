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
date: 12 January 2024
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

## Generating networks
To test phylogenetic network methods, one either needs to source or create a test set of networks. Creating them is often the simpler option, so methods to randomly generate phylogenetic networks are ready at hand. Moreover, these methods are often based on evolutionary models that are defined on a high level, i.e. with explicit events for  with processes such as speciation, extinction, and hybridization.

The paper [@janssen2021comparing] contains a comparison of several 'generators', including several previously existing ones (e.g., [@pons2019generation] and [@zhang2018bayesian]) and a new extention of a tree generator to networks.

# PhyloX Functionality



# Statement of Need
Currently, there is no Python package that enables a full workflow for analysing properties and methods of phylogenetic networks. Isolated scripts for this purpose do appear on GitHub or as pseudo-code regularly, most often as part of publications studying one method or one property [@zhang2018bayesian; @pons2019generation; @janssen2020combining; @janssen2020linear; @janssen2021rearranging]. Combining such scripts requires substantial work, for example because the phylogenetic networks themselves are represented by different Python classes with their own methods.

This package, PhyloX, aims to bring these scripts together: it standardizes implementations of several basic objects related to phylogenetic networks, such as the networks themselves, the labelling of the nodes, and rearrangement moves. It currently implements a limited but important set of basic functions: I/O for networks (e.g. lists of edges and extended newick format), network generation for test sets, comparing networks resulting from reconstruction methods, and computing several well-used network properties such as the reticulation number, the level, and the number of cherries.


# Related packages
As mentioned above, there are currently no Python packages that enable a complete workflow for phylogenetic networks. There are, however, a few Python packages that seem to fit that bill to a certain extent. We will argue that, despite these packages being available, there is still a need, or at least a great benefit, of using PhyloX.

## PhyloNetwork
Like PhyloX, [PhyloNetwork](https://github.com/bielcardona/PhyloNetwork) is a Python package based on NetworkX. It has a richer implementation for phylogenetic trees than PhyloX. For example, it includes more tree-specific rearrangement moves, the calculation of node properties such as the latest common ancestor (LCA), and some presets for drawing networks.

However, it has very few methods for phylogenetic networks and most of those methods are also included in PhyloX. Besides the network methods implemented in PhyloNetwork, PhyloX also includes isomorphism checking, rearrangement methods and distances for networks, more random network and tree generators, some functionality to combine trees or networks, and calculation of network properties. Another advantage of using PhyloX over PhyloNetwork is the inclusion of explicit random seeds. This is an important factor for the reproducibility of research.

Note that code from PhyloNetwork and PhyloX may be easy to combine, as both use NetworkX to implement the phylogenetic network class.

## Biopython - Phylo
This phylogenetics module, [Phylo](https://biopython.org/wiki/Phylo) [@talevich2012bio], of the Biopython package [@10.1093/bioinformatics/btp163] is built for phylogenetic analyses in Python. However, it is set up for phylogenetic trees only. The encoding of trees as sets of [clades](https://github.com/biopython/biopython/blob/17a9a5e41bafd9a85df18d4210e9293707d9e369/Bio/Phylo/PhyloXML.py#L321) does not easily allow extension to networks, which makes it unsuitable to use for these phylogenetic networks methods.

## DendroPy
Like Biopython's phylogenetics package, the [DendroPy](https://github.com/jeetsukumaran/DendroPy) package focuses on phylogenetic trees [@sukumaran2010dendropy]. Unlike Biopython, the implementation of the trees in DendroPy does seem to be graph based, making it more feasible to implement phylogenetic networks in DendroPy. This could still require large changes, as some properties of trees are built into the code on a fairly fundamental level, such as each node having (at most one) [parent node](https://github.com/jeetsukumaran/DendroPy/blob/cc82ab774ed83831b5c5125278d88c3c614c2d8a/src/dendropy/datamodel/treemodel/_node.py#L55C14-L55C26).


# Availability
The code of PhyloX is available as an open source project on [GitHub](https://github.com/RemieJanssen/PhyloX) under the BSD 3-Clause license. The package is also available via [PyPI](https://pypi.org/project/phylox/), so it can be installed via pip (or pip in conda), and updates to the release branch are automatically converted into new versions of the package. The releases are recorded in [Zenodo](https://zenodo.org/records/10122404) so persistent identifiers can be used to cite specific releases of the software. When citing this software, please make sure to also cite the original source of the code, which is mentioned in the [documentation](https://phylox.readthedocs.io/) of each method or class.

\section{Acknowledgements}
Most of the code has been written in the form of separate scripts during the author's PhD project, which was conducted under Leo van Iersel's \orcidlink{0000-0001-7142-4706} Vidi grant: 639.072.6

Anyone willing to contribute is very welcome to do so via pull requests and issues on github!

# References
