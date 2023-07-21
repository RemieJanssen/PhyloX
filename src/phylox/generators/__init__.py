"""
A module for generating (random) phylogenetic networks.

Several types of generators have been implemented, the algorithms are based on work by other authors.
These specific implementations for (heath, lgt, mcmc, trees, zods) are adapted from the ones used for the paper
"Comparing the topology of phylogenetic network generators" 
by Remie Janssen and Pengy Liu (2021).

The randomTC generator is based on the algorithm originally by Yukihiro Murakami, and improved by Remie Janssen to create test sets for the paper
"Linear Time Algorithm for Tree-Child network Containment" by Remie Janssen and Yukihiro Murakami (2020).

The beta-splitting tree generator is adapted from a script provided by a (colleague of) Pengyu Liu.
Which in turn was based on the Beta-splitting model by Aldous (1996).

The lgt generator is based on the algorithm from:
Joan Carles Pons, Celine Scornavacca, and Gabriel Cardona. Generation of level-k LGT
networks. IEEE/ACM transactions on computational biology and bioinformatics, 17(1):158–164,
2019.

The ZODS generator is based on the algorithm from:
Zhang, C., Ogilvie, H.A., Drummond, A.J., Stadler, T.: Bayesian inference of species networks from multilocus
sequence data. Molecular biology and evolution 35(2), 504–517 (2018)
"""
