# Supplementary files for "Generation of level-k LGT networks"
## By J.C. Pons, C. Scornavacca, and G. Cardona

The requirements to run these jupyter notebooks are python 3.6 with the library `phylonetwork` installed (can be installed via `pip`). The library `matplotlib` is also necessary to reproduce the plots in `evolution_model.py`

The included files are:

* `counting_generators.ipynb`: This script generates all the different LGT generators of (small) given level up to isomorphism. It outputs the number of such generators.
* `evolution_model.ipynb`: This script implements the stochastic model of evolution described in the paper. It performs several simulations with different parameters and displays the results.
* `Lev5Gen.txt`: This file contains the level-k generators (k=1,...,5) obtained from http://igm.univ-mlv.fr/~gambette/ProgGenerators.php in DOT format.
