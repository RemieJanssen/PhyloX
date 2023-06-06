def is_binary(network):
    binary_node_types = [
        [0, 1],  # root
        [0, 2],  # root
        [1, 2],  # tree node
        [2, 1],  # reticulation
        [1, 0],  # leaf
    ]
    for node in network.nodes:
        degrees = [network.in_degree(node), network.out_degree(node)]
        if degrees not in binary_node_types:
            return False
    return True

