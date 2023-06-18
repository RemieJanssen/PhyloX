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
    """

    new_node = -1
    while new_node in network.nodes or new_node in exclude:
        new_node -= 1
    return new_node
