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
 - name: Rijksinstituut voor Volksgezondheid en Milieu
   index: 1
date: 12 January 2024
bibliography: paper.bib
---

# Summary

PhyloX is a Python package with tools for generating, manipulating, and analyzing phylogenetic networks. It uses the NetworkX package [@SciPyProceedings_11] for basic graph operations. This has the side effect that the powerful graph tools from networkX can be used directly on the phylogenetic networks as well. The aim of the package is to be of general use to phylogenetic network researchers, but the focus is currently skewed towards random generation of networks, cherry-picking methods, rearrangement operations, and the identification of classes of networks.

# Statement of Need
Currently, there is no Python package that enables a full workflow for analysing properties and methods of phylogenetic networks. Isolated scripts for this purpose do appear on GitHub or as pseudo-code regularly, as part of publications studying one method or one property \citep[e.g.][]{zhang2018bayesian, pons2019generation, janssen2020combining, janssen2020linear, janssen2021rearranging}. Combining such scripts into one requires quite some work, for example because the phylogenetic networks themselves are represented by different Python classes with their own methods.

This package, PhyloX, aims to bring these scripts together: it standardizes implementations of several basic objects related to pylogenetic networks, such as the networks themselves, the labelling of the nodes, and rearrangement moves. It currently implements a limited but important set of basic functions: I/O for networks (e.g. lists of edges and extended newick format), network generation for test sets, comparing networks resulting from reconstruction methods, and computing several basic network properties.

# Related packages
In the statement of need we state that there are no pyton packages that enable a complete workflow for phylogenetic networks. There are, however, a few Python packages that seem to fit that bill. We will argue that, despite these packages being available, there is still a need or at least a great benefit of using PhyloX.

## PhyloNetwork
[PhyloNetwork](https://github.com/bielcardona/PhyloNetwork) is also a Python package based on NetworkX. PhyloNetwork has a richer implementation for phylogenetic trees than PhyloX. For example, it includes more tree-specific rearragement moves, the calculation of node properties such as the latest common ancestor (LCA), and some presets for drawing networks. 

However, it has very few methods for phylogenetic networks and most methods in PhyloNetwork are also included in PhyloX. Next to the network methods implemented in PhyloNetwork, PhyloX also includes isomorphism checking, rearrangement methods and distances for networks, more random network and tree generators, some functionality to combine trees or networks, and calculation of network properties. Another advantage of using PhyloX over PhyloNetwork is the inclusion of explicit random seeds. This is an important factor for the reproducibility of research.

Note that code from PhyloNetwork and PhyloX may be easy to combine, as both use NetworkX to implement the phylogenetic network class.

## Biopython - Phylo

This phylogenetics part, [Phylo](https://biopython.org/wiki/Phylo) [@talevich2012bio], of the Biopython package [@10.1093/bioinformatics/btp163] is built for phylogenetic analyses in Python. However, it set up for phylogenetic trees only. The encoding of trees as sets of [clades](https://github.com/biopython/biopython/blob/17a9a5e41bafd9a85df18d4210e9293707d9e369/Bio/Phylo/PhyloXML.py#L321) does not easily allow extension to networks, which makes it unsuitable to use for these phylogenetic networks methods.

## DendroPy
Like Biopython's phylogenetics package and as its name suggests, the [DendroPy](https://github.com/jeetsukumaran/DendroPy) package focuses on phylogenetic trees [@sukumaran2010dendropy]. Unlike Biopython, the implementation of the trees in DendroPy do seem to be more graph based, making it more feasible to implement phylogenetic networks in Dendropy. This may still require large changes, as some properties of trees are built into the code on a fairly fundamental level, such as each node having (at most one) [parent node](https://github.com/jeetsukumaran/DendroPy/blob/cc82ab774ed83831b5c5125278d88c3c614c2d8a/src/dendropy/datamodel/treemodel/_node.py#L55C14-L55C26). 

# Availability
The code of PhyloX is available as open source project on [GitHub](https://github.com/RemieJanssen/PhyloX) under the BSD-3-Clause license. The package is also available via [PyPI](https://pypi.org/project/phylox/), so it can be installed via pip (or pip in conda), and updates to the release branch are automatically converted into new versions of the package. The releases are recorded in [Zenodo](https://zenodo.org/records/10122404) so persistent identifiers can be used to cite specific releases of the software. When citing this software, please make sure to also cite the original source of the code, which is mentioned in the [documentation](https://phylox.readthedocs.io/) of each method or class.

\section{Acknowledgements}
Most of the code has been written in the form of separate scripts during the author's PhD project, which was conducted under Leo van Iersel's \orcidlink{0000-0001-7142-4706} Vidi grant: 639.072.6

Anyone willing to contribute is very welcome to do so via pull requests and issues on github!

# References