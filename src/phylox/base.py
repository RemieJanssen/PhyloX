from copy import deepcopy

from phylox.constants import LABEL_ATTR, LENGTH_ATTR


def find_unused_node(network, exclude=[]):
    """
    Find an unused node in the network.

    Parameters
    ----------
    network : networkx.DiGraph
        The network to find an unused node in.
    exclude : list
        A list of additional nodes to exclude from the search.

    Returns
    -------
    int
        The unused node.

    Examples
    --------
    >>> import networkx as nx
    >>> import phylox
    >>> network = nx.DiGraph()
    >>> network.add_edges_from([(0, 1), (1, 2), (2, 3)])
    >>> phylox.find_unused_node(network)
    -1
    >>> phylox.find_unused_node(network, exclude=[-1])
    -2
    """

    new_node = -1
    while new_node in network.nodes or new_node in exclude:
        new_node -= 1
    return new_node


def suppress_node(network, node):
    """
    Suppresses a degree-2 node node and returns true if successful.
    The new arc has length length(p,node)+length(node,c).
    Returns false if node is not a degree-2 node.

    Parameters
    ----------
    network : phylox.DiNetwork
        The network to suppress a node in.
    node : str or int
        The node to suppress.

    Returns
    -------
    bool
        True if the node was suppressed, False otherwise.

    Examples
    --------
    >>> from phylox import DiNetwork
    >>> network = DiNetwork()
    >>> network.add_edges_from([(0, 1, {'length': 1}), (1, 2, {'length': 1})])
    >>> suppress_node(network, 1)
    True
    >>> network.edges(data=True)
    OutEdgeDataView([(0, 2, {'length': 2})])
    """
    if not (network.out_degree(node) == 1 and network.in_degree(node) == 1):
        return False
    parent = network.parent(node)
    child = network.child(node)
    network.add_edges_from([(parent, child, network[parent][node])])
    if LENGTH_ATTR in network[parent][node] and LENGTH_ATTR in network[node][child]:
        network[parent][child][LENGTH_ATTR] = (
            network[parent][node][LENGTH_ATTR] + network[node][child][LENGTH_ATTR]
        )
    network.remove_node(node)
    return True


def remove_unlabeled_leaves(network, inplace=True):
    """
    Iteratively removes unlabeled leaves until none are left, then suppresses all degree-2 nodes.

    Parameters
    ----------
    network : phylox.DiNetwork
        The network to remove unlabeled leaves from.
    inplace : bool
        Whether to modify the network in place or return a copy.

    Returns
    -------
    phylox.DiNetwork
        The network with unlabeled leaves removed.

    Examples
    --------
    >>> from phylox import DiNetwork
    >>> from phylox.constants import LABEL_ATTR
    >>> network = DiNetwork(
    ...     edges=[(0, 1), (1, 2), (1, 3), (3, 4), (3, 5)],
    ...     labels=[(2, 'a'), (4, 'b')],
    ... )
    >>> network = remove_unlabeled_leaves(network)
    >>> network.edges()
    OutEdgeView([(0, 1), (1, 2), (1, 4)])
    >>> network.nodes(data=True)
    NodeDataView({0: {}, 1: {}, 2: {'label': 'a'}, 4: {'label': 'b'}})
    """
    if not inplace:
        network = deepcopy(network)

    nodes_to_check = set(deepcopy(network.nodes))
    while nodes_to_check:
        v = nodes_to_check.pop()
        if network.out_degree(v) == 0 and LABEL_ATTR not in network.nodes[v]:
            for p in network.predecessors(v):
                nodes_to_check.add(p)
            network.remove_node(v)
    list_nodes = deepcopy(network.nodes)
    for v in list_nodes:
        suppress_node(network, v)
    return network
