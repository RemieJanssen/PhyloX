"""
The module handles parsing of Newick strings.

The edges can have attached properties, each edge can have 0, 1, or 3 properties. 
If there is just one, it will be read as a branch length.
If there are three, then the first is the branch length, the second is the bootstrap value, and the third is the inheritance probability along that edge (useful for incoming edges of reticulation nodes).
"""

import json
import re
from copy import deepcopy

from phylox import DiNetwork
from phylox.constants import LABEL_ATTR, LENGTH_ATTR, PROBABILITY_ATTR, RETIC_PREFIX


def dinetwork_to_extended_newick(network, simple=False):
    """
    Converts a phylogenetic network to a Newick string.
    The newick string has :length:bootstrap:probability annotations if any edge has a bootstrap or probability.
    If only lengths are available, the newick string has :length annotations.

    :param network: a phylogenetic network, i.e., a phylox DiNetwork.
    :param simple: Boolean, indicating whether to create a simple newick string without parameters
    :return: a string in extended Newick format for phylogenetic networks.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.constants import LENGTH_ATTR
    >>> from phylox.newick_parser import dinetwork_to_extended_newick
    >>> network = DiNetwork(
    ...     edges=[
    ...         (0, 1, {LENGTH_ATTR: 1.0}),
    ...         (0, 2, {LENGTH_ATTR: 1.0}),
    ...         (1, 3, {LENGTH_ATTR: 1.0}),
    ...         (2, 3, {LENGTH_ATTR: 1.0}),
    ...         (1, 4, {LENGTH_ATTR: 1.0}),
    ...         (2, 5, {LENGTH_ATTR: 1.0}),
    ...         (3, 6, {LENGTH_ATTR: 1.0}),
    ...     ],
    ...     labels=((0, "A"), (4, "B"), (5, "C"), (6, "D"), (3, "E")),
    ... )
    >>> newick = dinetwork_to_extended_newick(network)
    >>> "C:1.0" in newick
    True
    """

    cut_network = deepcopy(network)

    roots = cut_network.roots
    if len(roots) > 1:
        raise ValueError("Network has more than one root.")
    root = roots.pop()

    for retic_id, node in enumerate(
        [node for node in cut_network.nodes if cut_network.is_reticulation(node)]
    ):
        _split_reticulation_node(cut_network, node, retic_id=retic_id)

    has_lengths = any(
        LENGTH_ATTR in cut_network[parent][child] for parent, child in cut_network.edges
    )
    has_bootstraps = any(
        "bootstrap" in cut_network[parent][child] for parent, child in cut_network.edges
    )
    has_probabilities = any(
        PROBABILITY_ATTR in cut_network[parent][child]
        for parent, child in cut_network.edges
    )

    def node_to_newick(node):
        node_label = cut_network.nodes[node].get(LABEL_ATTR, "")

        if cut_network.is_leaf(node):
            return node_label

        children_strings = []
        for child in cut_network.successors(node):
            child_str = str(node_to_newick(child))
            if (has_bootstraps or has_probabilities) and not simple:
                length = cut_network[node][child].get(LENGTH_ATTR, "")
                bootstrap = cut_network[node][child].get("bootstrap", "")
                probability = cut_network[node][child].get(PROBABILITY_ATTR, "")
                child_str += f":{length}:{bootstrap}:{probability}"
            elif has_lengths and not simple:
                child_str += f":{cut_network[node][child].get(LENGTH_ATTR, '')}"
            children_strings.append(child_str)

        return "(" + ",".join(children_strings) + ")" + node_label

    newick = node_to_newick(root)
    return newick + ";"


def _split_reticulation_node(network, node, retic_id):
    """
    Splits a reticulation node into multiple nodes.

    :param network: a phylogenetic network, i.e., a phylox DiNetwork.
    :param node: a node of network, indicating the reticulation node to be split.
    :return: a phylogenetic network, i.e., a phylox DiNetwork.

    :note: This function is used by dinetwork_to_extended_newick. It modifies the network in place.
    """

    parents = [
        (parent, network[parent][node].get(PROBABILITY_ATTR, 0))
        for parent in network.predecessors(node)
    ]
    parents.sort(key=lambda x: -x[1])
    keep_parent, keep_probability = parents[0]
    node_label = network.nodes[node].get(LABEL_ATTR, "")
    network.nodes[node][LABEL_ATTR] = node_label + "#R" + str(retic_id)
    new_node_label = node_label + "#H" + str(retic_id)
    for parent, probability in parents[1:]:
        new_node = network.find_unused_node()
        network.add_edge(parent, new_node)
        network.nodes[new_node][LABEL_ATTR] = new_node_label

        length = network[parent][node].get(LENGTH_ATTR)
        if length is not None:
            network[parent][new_node][LENGTH_ATTR] = length
        bootstrap = network[parent][node].get("bootstrap")
        if bootstrap is not None:
            network[parent][new_node]["bootstrap"] = bootstrap
        probability = network[parent][node].get(PROBABILITY_ATTR)
        if probability is not None:
            network[parent][new_node][PROBABILITY_ATTR] = probability
        network.remove_edge(parent, node)
    return network


def extended_newick_to_dinetwork(newick, internal_labels=False):
    """
    Converts a Newick string to a networkx DAG with leaf labels.
    The newick string may or may not have length:bootstrap:probability annotations.
    The newick string may or may not have internal node labels.
    The newick string may or may not have hybrid nodes.

    :param newick: a string in extended Newick format for phylogenetic networks.
    :param internal_labels: a boolean, indicating whether the internal nodes of the network are labeled.
    :return: a phylogenetic network, i.e., a networkx digraph with leaf labels represented by the `label' node attribute.

    :example:
    >>> newick = "(A:1.1,B:1.2,(C:1.3,D:1.4)E:1.6)F;"
    >>> network = extended_newick_to_dinetwork(newick)
    >>> {network.nodes[leaf].get("label") for leaf in network.leaves} == {'A', 'B', 'C', 'D'}
    True
    >>> node_for_label_A = network.label_to_node_dict['A']
    >>> p = network.parent(node_for_label_A)
    >>> network[p][node_for_label_A]['length']
    1.1
    """

    network = DiNetwork()
    network_json = _newick_to_json(newick)[0]
    network = _json_to_dinetwork(network_json, network=network)
    return network


def _newick_to_json(newick):
    """
    Converts a newick string to a json representing the nodes in a nested manner.

    :param newick: a string in newick format for phylogenetic networks.
    :return: a json representing the nodes in a nested manner.

    :note: This function is used by extended_newick_to_dinetwork. It modifies the network in place.
    """

    nested_list = [{"children": [], "label_and_attr": ""}]
    while newick:
        character = newick[0]
        newick = newick[1:]
        if character == "(":
            newick, child = _newick_to_json(newick)
            nested_list[-1]["children"] += child
        elif character == ")":
            return newick, nested_list
        elif character == ",":
            nested_list += [{"children": [], "label_and_attr": ""}]
        elif character == ";":
            pass
        else:
            nested_list[-1]["label_and_attr"] += character
    return nested_list


def _json_to_dinetwork(json, network=None, root_node=None):
    """
    Converts a json string to a phylox DiNetwork

    :param newick_json: a string in json format for phylogenetic networks.
    :param network: a phylogenetic network, i.e., a phylox DiNetwork.
    :param root_node: a node of network, indicating the root of the network.
    :return: a phylogenetic network, i.e., a phylox DiNetwork.

    :note: This function is used by extended_newick_to_dinetwork. It modifies the network in place.
    """
    network = network or DiNetwork()

    node_attrs = _label_and_attrs_to_dict(json["label_and_attr"])
    node = json.get("retic_id") or root_node or network.find_unused_node()
    network.add_node(node)
    if LABEL_ATTR in node_attrs:
        network.nodes[node][LABEL_ATTR] = node_attrs[LABEL_ATTR]
    for child_dict in json.get("children", []):
        child_attrs = _label_and_attrs_to_dict(child_dict["label_and_attr"])
        child_attrs_without_label_and_children = {
            k: v
            for k, v in child_attrs.items()
            if k not in (LABEL_ATTR, "children", "retic_id")
        }
        child = child_attrs.get("retic_id", network.find_unused_node())
        network.add_edge(node, child, **child_attrs_without_label_and_children)
        _json_to_dinetwork(child_dict, network, root_node=child)
    return network


def _label_and_attrs_to_dict(label_and_attrs):
    """
    converts the label and attr part of an extended newick string
    for one node to a dictionary.
    For example, the string "A:1.1:0.9:0.8" is converted to
    {"label": "A", "length": 1.1, "bootstrap": 0.9, "probability": 0.8}
    """
    attrs_dict = {LABEL_ATTR: label_and_attrs}
    if ":" in label_and_attrs:
        label = label_and_attrs.split(":")[0]
        attrs = label_and_attrs.split(":")[1:]
        if len(attrs) == 1:
            attrs_dict = {
                LABEL_ATTR: label,
                LENGTH_ATTR: float(attrs[0]),
            }
        elif len(attrs) == 3:
            attrs_dict = {LABEL_ATTR: label}
            if attrs[0]:
                attrs_dict[LENGTH_ATTR] = float(attrs[0])
            if attrs[1]:
                attrs_dict["bootstrap"] = float(attrs[1])
            if attrs[2]:
                attrs_dict[PROBABILITY_ATTR] = float(attrs[2])

    if "#" in attrs_dict[LABEL_ATTR]:
        label, retic_id = attrs_dict[LABEL_ATTR].split("#")
        attrs_dict[LABEL_ATTR] = label
        attrs_dict["retic_id"] = RETIC_PREFIX + retic_id[1:]
    if LABEL_ATTR in attrs_dict and attrs_dict[LABEL_ATTR] == "":
        attrs_dict.pop(LABEL_ATTR)
    return attrs_dict
