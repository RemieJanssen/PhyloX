from copy import deepcopy
from enum import Enum

from phylox.base import find_unused_node, suppress_node
from phylox.constants import LABEL_ATTR, LENGTH_ATTR


class CHERRYTYPE(Enum):
    CHERRY = 1
    RETICULATEDCHERRY = 2
    NONE = 0


def find_all_reducible_pairs(network):
    """
    Finds all reducible pairs (cherries and reticulated cherries) in the
    network.

    :param network: a phylogenetic network.
    :return: a set of reducible pairs (cherries and reticulated cherries) in the network.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.cherrypicking.base import find_all_reducible_pairs
    >>> network = DiNetwork(
    ...     edges=[(-1,0),(0,1),(1,2),(1,3),(2,3),(2,4),(3,5),(0,6),(6,7),(6,8)],
    ... )
    >>> reducible_pairs = find_all_reducible_pairs(network)
    >>> reducible_pairs == {(7,8),(8,7),(5,4)}
    True
    """
    reducible_pairs = set()
    for l in network.leaves:
        reducible_pairs = reducible_pairs.union(
            find_reducible_pairs_with_second(network, l)
        )
    return reducible_pairs


def find_reducible_pairs_with_second(N, x):
    """
    Finds a list of reducible pairs (cherries and reticulated cherries) in the
    network N with leaf x as second element of the pair.

    :param N: a phylogenetic network.
    :param x: a leaf of the network N.
    :return: a list of reducible pairs (cherries and reticulated cherries) in the network N with leaf x as second element of the pair.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.cherrypicking.base import find_reducible_pairs_with_second
    >>> network = DiNetwork(
    ...     edges=[(-1,0), (0,1), (0,2), (1,2), (1,3), (2,4)],
    ... )
    >>> find_reducible_pairs_with_second(network, 3)
    [(4, 3)]
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
            sibling_child = N.child(sibling)
            if N.out_degree(sibling_child) == 0:
                reducible_pairs.append((sibling_child, x))
    return reducible_pairs


def find_reducible_pairs_with_first(N, x):
    """
    Finds a list of reducible pairs (cherries and reticulated cherries) in the
    network N with leaf x as first element of the pair.

    :param N: a phylogenetic network.
    :param x: a leaf of the network N.
    :return: a list of reducible pairs (cherries and reticulated cherries) in the network N with leaf x as first element of the pair.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.cherrypicking.base import find_reducible_pairs_with_first
    >>> network = DiNetwork(
    ...     edges=[(-1,0), (0,1), (0,2), (1,2), (1,3), (2,4)],
    ... )
    >>> find_reducible_pairs_with_first(network, 4)
    [(4, 3)]
    """
    if not N.is_leaf(x):
        raise ValueError("x must be a leaf of N")

    parent = N.parent(x)

    if N.is_tree_node(parent):
        return find_cherries_with_first(N, x)
    if N.is_reticulation(parent):
        return find_reticulated_cherry_with_first(N, x)
    else:
        return []


def find_reticulated_cherry_with_first(N, x):
    """
    Finds a list of reticulated cherries in the network N with leaf x as first
    element of the pair.

    Parameters
    ----------
    N : phylox.DiNetwork
        The network to find reticulated cherries in.
    x : str or int
        The leaf to find reticulated cherries with.

    Returns
    -------
    list
        A list of reticulated cherries in the network with leaf x as first
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


def find_cherries_with_first(network, x):
    """
    Finds a set of cherries in the network N with leaf x as first element of
    the pair.

    Parameters
    ----------
    network : phylox.DiNetwork
        The network to find cherries in.
    x : str or int
        The leaf to find cherries with.

    Returns
    -------
    set
        A set of cherries in the network with leaf x as first element of the
        pair.
    """
    cherries = set()
    parent = network.parent(x)
    for sibling in network.successors(parent):
        if sibling in network.leaves and sibling != x:
            cherries.add((x, sibling))
    return cherries


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


def reduce_pair(network, x, y, inplace=False, nodes_by_label=False):
    """
    Reduces the reducible pair (x,y) in the network.
    Note: Cache of network properties is not updated.

    Parameters
    ----------
    network : phylox.DiNetwork
        The network to reduce the reducible pair in.
    x : str or int
        The first element of the reducible pair.
    y : str or int
        The second element of the reducible pair.
    inplace : bool
        If True, the network is modified in place.
    nodes_by_label : bool
        If True, the nodes x and y are interpreted as labels.

    Returns
    -------
    phylox.DiNetwork
        The network with the reducible pair reduced.
    CHERRYTYPE
        The type of the reducible pair.

    Raises
    ------
    ValueError
        If x or y are not in the network.

    Examples
    --------
    >>> from phylox import DiNetwork
    >>> from phylox.cherrypicking.base import reduce_pair, CHERRYTYPE
    >>> network = DiNetwork(
    ...     edges=[(-1,0), (0,1), (0,2), (1,2), (1,3), (2,4)],
    ... )
    >>> network, cherry_type = reduce_pair(network, 4, 3)
    >>> cherry_type == CHERRYTYPE.RETICULATEDCHERRY
    True
    >>> set(network.edges) == {(-1, 0), (0, 3), (0, 4)}
    True
    """

    if not inplace:
        network = deepcopy(network)
    if nodes_by_label:
        x = network.labels[x][0]
        y = network.labels[y][0]

    cherry_type = check_reducible_pair(network, x, y)
    if cherry_type == CHERRYTYPE.CHERRY:
        px = network.parent(x)
        network.remove_node(x)
        suppress_node(network, px)
    if cherry_type == CHERRYTYPE.RETICULATEDCHERRY:
        px = network.parent(x)
        py = network.parent(y)
        network.remove_edge(py, px)
        suppress_node(network, px)
        suppress_node(network, py)
    return network, cherry_type


def check_reducible_pair(network, x, y):
    """
    Checks whether the pair (x,y) is a reducible pair in the network.

    :param network: a phylogenetic network.
    :param x: a leaf of the network.
    :param y: a leaf of the network.
    :return: the type of reducible pair (x,y) in the network.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.cherrypicking.base import check_reducible_pair, CHERRYTYPE
    >>> network = DiNetwork(
    ...     edges=[(-1,0), (0,1), (0,2), (1,2), (1,3), (2,4)],
    ... )
    >>> check_reducible_pair(network, 4, 3) == CHERRYTYPE.RETICULATEDCHERRY
    True
    """
    if not network.has_node(x):
        return CHERRYTYPE.NONE
    if not network.has_node(y):
        return CHERRYTYPE.NONE
    px = network.parent(x)
    py = network.parent(y)
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


# TODO: make work with cps with labels instead of node indices
def get_indices_of_reducing_pairs(sequence, network):
    """
    Checks which pairs of a sequence actually reduce a given network
    for a given cherry-picking sequence `sequence' reduces a given tree `tree'
    input:
        sequence: a list of pairs of leaves
        network: a network
    output:
        if the network is reduced by the sequence, returns the list of all indices of pairs that reduce the network
        otherwise returns False
    """
    network_copy = deepcopy(network)
    indices = []
    for i, pair in enumerate(sequence):
        network_copy, cherry_type = reduce_pair(network_copy, *pair)
        if cherry_type != CHERRYTYPE.NONE:
            indices += [i]
        if len(network_copy.edges) <= 1:
            return indices
    return False


def add_roots_to_sequence(sequence, reduced_trees_per_pair):
    """
    Modifies a cherry-picking sequence so that it represents a network with exactly one root.
    A sequence may be such that reconstructing a network from the sequence results in multiple roots
    This function adds some pairs to the sequence so that the network has a single root.
    args:
        sequence: the sequence to modify
        reduced_trees_per_pair: the sets of trees reduced by each pair in the sequence
    returns:
        the new sequence, and also the sets of trees reduced by each pair in the sequence, modified so that the new pairs are also represented (they reduce no trees)
    """
    leaves_encountered = set()
    roots = set()
    # The roots can be found by going back through the sequence and finding pairs where the second element has not been encountered in the sequence yet
    for pair in reversed(sequence):
        if pair[1] not in leaves_encountered:
            roots.add(pair[1])
        leaves_encountered.add(pair[0])
        leaves_encountered.add(pair[1])
    i = 0
    roots = list(roots)
    # Now add some pairs to make sure each second element is already part of some pair in the sequence read backwards, except for the last pair in the sequence
    for i in range(len(roots) - 1):
        sequence.append((roots[i], roots[i + 1]))
        # none of the trees are reduced by the new pairs.
        reduced_trees_per_pair.append(set())
        i += 1
    return sequence, reduced_trees_per_pair


def has_cherry(network, x, y):
    """
    Checks whether the pair (x,y) forms a cherry in the network

    Parameters
    ----------
    network : phylox.DiNetwork
        The network in which we want to check whether (x,y) is a cherry
    x : string
        The first element of the pair
    y : string
        The second element of the pair

    Returns
    -------
    bool
        True if (x,y) is a cherry in the network, False otherwise
    """
    if (not x in network.leaves) or (not y in network.leaves):
        return False
    px = network.parent(x)
    py = network.parent(y)
    return px == py


def cherry_height(network, x, y):
    """
    Returns the height of (x,y) if it is a cherry:
        i.e.: length(p,x)+length(p,y)/2
    Returns false otherwise

    Parameters
    ----------
    network : phylox.DiNetwork
        The network in which we want to check the height of cherry (x,y)
    x : string
        The first element of the pair
    y : string
        The second element of the pair

    Returns
    -------
    float
        The height of the cherry (x,y) if it is a cherry, False otherwise
    """
    print("cherry_height")
    print(network.edges(data=True))

    if (not x in network.leaves) or (not y in network.leaves):
        return False
    px = network.parent(x)
    py = network.parent(y)
    if px == py:
        height = [network[px][x][LENGTH_ATTR], network[py][y][LENGTH_ATTR]]
        return height
    if (py, px) in network.edges:
        height = [
            network[px][x][LENGTH_ATTR] + network[py][px][LENGTH_ATTR],
            network[py][y][LENGTH_ATTR],
        ]
        return height
    raise ValueError("x and y are not in the same cherry")


class CherryPickingMixin:
    @classmethod
    def from_cherry_picking_sequence(cls, sequence, heights=None, label_leaves=True):
        network = cls()
        heights = heights or [[1, 1]] * len(sequence)
        for pair, height in zip(reversed(sequence), reversed(heights)):
            add_pair(
                network, *pair, height=height, inplace=True, nodes_by_label=label_leaves
            )
        network._clear_cached()
        return network
