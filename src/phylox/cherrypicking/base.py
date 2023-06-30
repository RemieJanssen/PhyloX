from copy import deepcopy
from enum import Enum

from phylox.base import find_unused_node
from phylox.constants import LABEL_ATTR, LENGTH_ATTR


class CHERRYTYPE(Enum):
    CHERRY = 1
    RETICULATEDCHERRY = 2
    NONE = 0


def find_reducible_pairs_with_second(N, x):
    """
    Finds a list of reducible pairs (cherries and reticulated cherries) in the
    network N with leaf x as second element of the pair.
    """
    if not N.is_leaf(x):
        raise ValueError("x must be a leaf of N")

    parent = N.parent(x)
    if N.in_degree(parent) == 0:
        return []

    reducible_pairs = list()
    for sibling in N.successors(parent):
        if sibling == x:
            continue
        sibling_out_degree = N.out_degree(sibling)
        if sibling_out_degree == 0:
            reducible_pairs.append((sibling, x))
        if sibling_out_degree == 1:
            for sibling_child in N.successors(sibling):
                if N.out_degree(sibling_child) == 0:
                    reducible_pairs.append((sibling_child, x))
    return reducible_pairs


def find_reticulated_cherry_with_first(N, x):
    """
    Finds a list of reticulated cherries in the network N with leaf x as first
    element of the pair.
    """

    if not N.is_leaf(x):
        raise ValueError("x must be a leaf of N")

    parent = N.parent(x)
    if parent is None:
        return []
    if not N.is_reticulation(parent):
        return []

    reticulated_cherries = list()
    for pp in N.predecessors(parent):
        for ppc in N.successors(pp):
            if ppc == parent or not N.is_leaf(ppc):
                continue
            reticulated_cherries.append((x, ppc))
    return reticulated_cherries


def is_second_in_reducible_pair(network, x):
    for node in network.predecessors(x):
        px = node
    for cpx in network.successors(px):
        if cpx != x:
            if network.out_degree(cpx) == 0:
                return (cpx, x)
            if network.out_degree(cpx) == 1:
                for ccpx in network.successors(cpx):
                    if network.out_degree(ccpx) == 0:
                        return (ccpx, x)
    return False


def reduce_pair(network, x, y, inplace=False):
    if not inplace:
        network = deepcopy(network)

    cherry_type = check_reducible_pair(network, x, y)
    if cherry_type == CHERRYTYPE.CHERRY:
        px = network.parent(x)
        network.remove_node(x)
        if network.out_degree(px) == 1:
            ppx = network.parent(px)
            network.remove_node(px)
            network.add_edge(ppx, y)
    if cherry_type == CHERRYTYPE.RETICULATEDCHERRY:
        px = network.parent(x)
        py = network.parent(y)
        network.remove_edge(py, px)
        if network.in_degree(px) == 1:
            ppx = network.parent(px)
            network.add_edge(ppx, x)
            network.remove_node(px)
        if network.out_degree(py) == 1:
            ppy = network.parent(py)
            network.add_edge(ppy, y)
            network.remove_node(py)
    if inplace:
        # TODO empty cache for network properties?
        pass
    return network, cherry_type


def check_reducible_pair(network, x, y):
    if network.has_node(x):
        if network.has_node(y):
            for px in network.predecessors(x):
                for py in network.predecessors(y):
                    if px == py:
                        return CHERRYTYPE.CHERRY
                    if network.out_degree(px) == 1:
                        if px in network.successors(py):
                            return CHERRYTYPE.RETICULATEDCHERRY
    return CHERRYTYPE.NONE


def add_pair(network, x, y, height=[1, 1], inplace=False, nodes_by_label=False):
    """
    Adds a pair to the network, using the construction from a cherry-picking sequence
    :param x: first element of the pair
    :param y: second element of the pair
    :param height: height of the pair
    :param inplace: if true, the network is modified in place, otherwise a copy is returned
    :param nodes_by_label: if true, the nodes are indexed by their label, otherwise by their index
    :return: the network with the pair added
    """
    print("add pair", x, y, height)
    if not inplace:
        network = deepcopy(network)

    # if the network is empty, create a cherry (x,y)
    if len(network.leaves) == 0:

        node_x = 2 if nodes_by_label else x
        node_y = 3 if nodes_by_label else y
        root = find_unused_node(network, exclude=[node_x, node_y])
        parent = find_unused_node(network, exclude=[root, node_x, node_y])
        network.add_weighted_edges_from(
            [
                (root, parent, 0),
                (parent, node_x, height[0]),
                (parent, node_y, height[1]),
            ],
            weight=LENGTH_ATTR,
        )
        if nodes_by_label:
            network.label_to_node_dict[x] = node_x
            network.label_to_node_dict[y] = node_y
            network.nodes[node_x][LABEL_ATTR] = x
            network.nodes[node_y][LABEL_ATTR] = y
        network._clear_cached()
        return network

    node_y = network.label_to_node_dict.get(y) if nodes_by_label else y
    # if y is not in the network raise an error, as there is no way to add the pair and get a phylogenetic network
    if node_y is None or node_y not in network.leaves:
        raise ValueError("y is not in the network")
    # else, add the pair to the existing network
    # get edge data for edges around y
    parent_node_y = network.parent(node_y)
    length_incoming_y = network[parent_node_y][node_y].get(LENGTH_ATTR)
    # no_of_trees_incoming_y = network[parent_node_y][node_y].get("no_of_trees")
    edge_data = dict()
    height_goal_x = height[0]
    if length_incoming_y is not None:
        if height[1] < length_incoming_y:
            height_pair_y_real = height[1]
        else:
            height_pair_y_real = length_incoming_y
            height_goal_x += height[1] - height_pair_y_real
        edge_data[LENGTH_ATTR] = height_pair_y_real
    # if no_of_trees_incoming_y is not None:
    #     edge_data["no_of_trees"] = no_of_trees_incoming_y + len(red_trees - current_trees)

    old_edge_data = network.edges[parent_node_y, node_y]
    old_edge_data[LENGTH_ATTR] = length_incoming_y - height_pair_y_real

    # add all edges around y
    network.remove_edge(parent_node_y, node_y)
    new_parent_of_y = find_unused_node(network)
    network.add_edges_from(
        [
            (parent_node_y, new_parent_of_y, old_edge_data),
            (new_parent_of_y, node_y, edge_data),
        ]
    )

    # Now also add edges around x
    node_x = (
        network.label_to_node_dict.get(x, find_unused_node(network))
        if nodes_by_label
        else x
    )
    # x is not yet in the network, so make a cherry (x,y)
    if node_x not in network.leaves:
        network.add_edge(
            new_parent_of_y,
            node_x,
            # no_of_trees=len(red_trees),
            length=height_goal_x,
        )
        if nodes_by_label:
            network.nodes[node_x][LABEL_ATTR] = x
            network.label_to_node_dict[x] = node_x
        network._clear_cached()
        return network

    # x is already in the network, so create a reticulate cherry (x,y)
    parent_node_x = network.parent(node_x)
    length_incoming_x = network[parent_node_x][node_x][LENGTH_ATTR]
    # no_of_trees_incoming_x = network[parent_node_x][node_x]["no_of_trees"]
    # if x is below a reticulation, and the height of the new pair is above the height of this reticulation, add the new hybrid arc to the existing reticulation
    if network.in_degree(parent_node_x) > 1 and length_incoming_x <= height_goal_x:
        network.add_edge(
            new_parent_of_y,
            parent_node_x,
            # no_of_trees=len(red_trees),
            length=height_goal_x - length_incoming_x,
        )
        # network[parent_node_x][node_x]["no_of_trees"] += len(red_trees)
        network._clear_cached()
        return network

    # create a new reticulation vertex above x to attach the hybrid arc to
    height_pair_x = min(height_goal_x, length_incoming_x)
    new_parent_of_x = find_unused_node(network)
    old_edge_data = network.edges[parent_node_x, node_x]
    old_edge_data[LENGTH_ATTR] = length_incoming_x - height_pair_x

    network.remove_edge(parent_node_x, node_x)
    network.add_edges_from(
        [
            (parent_node_x, new_parent_of_x, old_edge_data),
            (
                new_parent_of_x,
                node_x,
                {LENGTH_ATTR: length_incoming_x - height_pair_x},
            ),  # "no_of_trees": no_of_trees_incoming_x + len(red_trees)
            (
                new_parent_of_y,
                new_parent_of_x,
                {LENGTH_ATTR: height_goal_x - height_pair_x},
            ),  # "no_of_trees": len(red_trees)
        ]
    )
    network._clear_cached()
    return network


class CherryPickingMixin:
    @classmethod
    def from_cherry_picking_sequence(cls, sequence, heights=None, label_leaves=True):
        network = cls()
        heights = heights or [[1, 1]] * len(sequence)
        for pair, height in zip(reversed(sequence), reversed(heights)):
            print("edges", network.edges(data=True))
            print("nodes", network.nodes(data=True))
            add_pair(
                network, *pair, height=height, inplace=True, nodes_by_label=label_leaves
            )
        network._clear_cached()
        print(sequence)
        print(network.edges)
        return network
