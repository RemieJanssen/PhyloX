"""
Adapted from script provided by (a colleage of) Pengyu Liu
which contains the code for Beta-splitting model (Aldous1996)

The beta-splitting model is a model for generating random binary trees.
The model is parameterized by a parameter beta > 0 which determines the shape of the tree.
"""

import random

import numpy as np
from scipy.special import loggamma

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR

############################################
# Simulation functions
############################################

# _a_n is a normalizing constant defined in
# Equation (2) of Aldous1996 (so the sum of
# the values is equal to 1. It is not
# computed here to save time.
def _a_n(beta):
    return 1


# Compute the "probability" to split n in (i|n-1), where i=1,..,n-1
def _compute_split_probability(n, beta):
    q_n = []
    for i in range(1, n):
        q_i_n = np.exp(
            (loggamma(beta + i + 1) + loggamma(beta + n - i + 1))
            - ((_a_n(beta) + loggamma(i + 1) + loggamma(n - i + 1)))
        )
        q_n.append(q_i_n)
    return q_n


# n: number of tips
def simulate_beta_splitting(n, beta):
    """
    Simulate a random binary tree with n leaves using the beta-splitting model.

    :param n: the number of leaves of the tree.
    :param beta: the beta parameter of the beta-splitting model.
    :return: a random binary tree with n leaves using the beta-splitting model.

    :example:
    >>> from phylox.generators.trees.beta_splitting_tree import simulate_beta_splitting
    >>> tree = simulate_beta_splitting(5, 1)
    >>> len(tree.leaves)
    5
    """
    # Initialize tree.
    tree = DiNetwork()
    tree.add_edge(-1, n + 1)
    tree.nodes[n + 1][LABEL_ATTR] = n
    last_internal_node = n + 1
    last_leaf_node = 0
    queue = [n + 1]
    # Insert one node at each iteration.
    while queue:
        node = queue.pop()
        n_node = tree.nodes[node].get(LABEL_ATTR)
        # Internal node. Splits again.
        if n_node > 1:
            # Compute the "probability" to split n in (i|n-1), where i=1,..,n-1
            q_n = _compute_split_probability(n_node, beta)
            split = random.choices(population=list(range(1, n_node)), weights=q_n, k=1)[
                0
            ]
            # Create children.
            for new_n in [split, n_node - split]:
                if new_n == 1:
                    tree.add_edge(node, last_leaf_node + 1)
                    last_leaf_node += 1
                else:
                    tree.add_edge(node, last_internal_node + 1)
                    tree.nodes[last_internal_node + 1][LABEL_ATTR] = new_n
                    queue.append(last_internal_node + 1)
                    last_internal_node += 1
    # Return tree.
    return tree
