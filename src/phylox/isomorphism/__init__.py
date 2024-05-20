"""
This module handles isomorphisms/automorphism of phylogenetic networks.

The functions in this module are used to check whether two networks are isomorphic (with or without labels).
In addition, it can count the number of automorphisms of a network, which is used in the mcmc network generator to correct for symmetries.
"""

from .base import is_isomorphic, count_automorphisms
