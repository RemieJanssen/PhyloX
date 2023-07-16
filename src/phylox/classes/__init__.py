"""
Classes for phylogenetic networks.

This module contains methods for checking whether a network is in a certain class.
Such as binary, tree-based, tree-child, orchard, or stack-free.

Originally written by Remie Janssen for the paper
"Comparing the topology of phylogenetic network generators"
by Remie Janssen and Pengy Liu.
"""


from .dinetwork import (
    is_binary,
    is_tree_based,
    is_tree_child,
    is_orchard,
    is_stack_free,
)
from .networkclass import *
