import random

import networkx as nx

from phylox.cherrypicking.base import CherryPickingMixin
from phylox.constants import LABEL_ATTR, LENGTH_ATTR


class DiNetwork(nx.DiGraph, CherryPickingMixin):
    def __init__(self, *args, **kwargs):
        edges = kwargs.get("edges", [])
        super().__init__(edges, *args, **kwargs)
        self.add_nodes_from(kwargs.get("nodes", []))
        self.label_to_node_dict = {}
        for label in kwargs.get("labels", []):
            self.nodes[label[0]][LABEL_ATTR] = label[1]
            self.label_to_node_dict[label[1]] = label[0]

    def _clear_cached(self):
        for attr in ["_leaves", "_reticulations", "_roots", "_reticulation_number, _labels"]:
            if hasattr(self, attr):
                delattr(self, attr)

    @classmethod
    def from_newick(cls, newick):
        pass

    def _set_leaves(self):
        self._leaves = set([node for node in self.nodes if self.is_leaf(node)])
        return self._leaves

    @property
    def leaves(self):
        if not hasattr(self, "_leaves"):
            self._set_leaves()
        return self._leaves

    def _set_reticulations(self):
        self._reticulations = set(
            [node for node in self.nodes if self.is_reticulation(node)]
        )
        return self._reticulations

    @property
    def reticulations(self):
        if not hasattr(self, "_retculations"):
            self._set_reticulations()
        return self._reticulations

    def _set_roots(self):
        self._roots = set([node for node in self.nodes if self.is_root(node)])
        return self._roots

    @property
    def roots(self):
        if not hasattr(self, "_roots"):
            self._set_roots()
        return self._roots

    @property
    def reticulation_number(self):
        if not hasattr(self, "_reticulation_number"):
            self._reticulation_number = sum(
                [max(self.in_degree(node) - 1, 0) for node in self.nodes]
            )
        return self._reticulation_number

    def _set_labels(self):
        self._labels = {}
        for node in self.nodes:
            if LABEL_ATTR in self.nodes[node]:
                label = self.nodes[node][LABEL_ATTR]
                if label not in self._labels:
                    self._labels[label] = []
                self._labels[self.nodes[node][LABEL_ATTR]] += [node]
        return self._labels

    @property
    def labels(self):
        if not hasattr(self, "_labels"):
            self._set_labels()
        return self._labels

    def child(self, node, exclude=[], randomNodes=False):
        """
        Finds a child node of a node.

        :param node: a node of self.
        :param exclude: a set of nodes of self.
        :param randomNodes: a boolean value.
        :return: a child of node that is not in the set of nodes exclude. If randomNodes, then this child node is selected uniformly at random from all candidates.
        """
        child = None
        for c in self.successors(node):
            if c not in exclude:
                if not randomNodes:
                    return c
                elif child is None or random.getrandbits(1):
                    # As there are at most two children, we can simply replace the previous child with probability .5 to get a random parent
                    child = c
        return child

    def parent(self, node, exclude=[], randomNodes=False):
        """
        Finds a parent of a node in a network.

        :param node: a node in the network.
        :param exclude: a set of nodes of the network.
        :param randomNodes: a boolean value.
        :return: a parent of node that is not in the set of nodes exclude. If randomNodes, then this parent is selected uniformly at random from all candidates.
        """
        parent = None
        for p in self.predecessors(node):
            if p not in exclude:
                if not randomNodes:
                    return p
                elif parent is None or random.getrandbits(1):
                    # As there are at most two parents, we can simply replace the previous parent with probability .5 to get a random parent
                    parent = p
        return parent

    def is_reticulation(self, node):
        return self.out_degree(node) <= 1 and self.in_degree(node) > 1

    def is_leaf(self, node):
        return self.out_degree(node) == 0 and self.in_degree(node) > 0

    def is_root(self, node):
        return self.in_degree(node) == 0

    def is_tree_node(self, node):
        return self.out_degree(node) > 1 and self.in_degree(node) <= 1
