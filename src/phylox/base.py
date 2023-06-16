def find_unused_node(network, exclude=[]):
    new_node = -1
    while new_node in network.nodes or new_node in exclude:
        new_node -= 1
    return new_node
