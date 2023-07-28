import re
import json
from phylox import DiNetwork
from phylox.base import find_unused_node
from phylox.constants import LABEL_ATTR, LENGTH_ATTR, RETIC_PREFIX

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
    network_json = newick_to_json(newick)[0]
    network = json_to_dinetwork(network_json, network=network)
    return network


def newick_to_json(newick):
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
            newick, child = newick_to_json(newick)
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


def json_to_dinetwork(json, network=None, root_node=None):
    """
    Converts a json string to a phylox DiNetwork

    :param newick_json: a string in json format for phylogenetic networks.
    :param network: a phylogenetic network, i.e., a phylox DiNetwork.
    :param root_node: a node of network, indicating the root of the network.
    :return: a phylogenetic network, i.e., a phylox DiNetwork.

    :note: This function is used by extended_newick_to_dinetwork. It modifies the network in place.
    """
    network = network or DiNetwork()

    node_attrs = label_and_attrs_to_dict(json["label_and_attr"])
    node = json.get("retic_id") or root_node or find_unused_node(network)
    network.add_node(node)
    if "label" in node_attrs:
        network.nodes[node]["label"] = node_attrs["label"]
    for child_dict in json.get("children", []):
        child_attrs = label_and_attrs_to_dict(child_dict["label_and_attr"])
        child_attrs_without_label_and_children = {
            k: v
            for k, v in child_attrs.items()
            if k not in ("label", "children", "retic_id")
        }
        child = child_attrs.get("retic_id", find_unused_node(network))
        network.add_edge(node, child, **child_attrs_without_label_and_children)
        json_to_dinetwork(child_dict, network, root_node=child)
    return network


def label_and_attrs_to_dict(label_and_attrs):
    """
    converts the label and attr part of an extended newick string
    for one node to a dictionary.
    For example, the string "A:1.1:0.9:0.8" is converted to
    {"label": "A", "length": 1.1, "bootstrap": 0.9, "probability": 0.8}
    """
    attrs_dict = {"label": label_and_attrs}
    if ":" in label_and_attrs:
        label = label_and_attrs.split(":")[0]
        attrs = label_and_attrs.split(":")[1:]
        if len(attrs) == 1:
            attrs_dict = {
                "label": label,
                "length": float(attrs[0]),
            }
        elif len(attrs) == 3:
            attrs_dict = {
                "label": label,
                "length": float(attrs[0]),
                "bootstrap": float(attrs[1]),
                "probability": float(attrs[2]),
            }
    if "#" in attrs_dict["label"]:
        label, retic_id = attrs_dict["label"].split("#")
        attrs_dict["label"] = label
        attrs_dict["retic_id"] = RETIC_PREFIX + retic_id[1:]
    if "label" in attrs_dict and attrs_dict["label"] == "":
        attrs_dict.pop("label")
    return attrs_dict
