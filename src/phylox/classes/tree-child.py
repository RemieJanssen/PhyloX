def is_tree_child(network):
    for node in network.nodes:
        if network.is_leaf(node):
            continue
        if all(
            [
                network.is_reticulation(child)
                for child in network.successors(node)
            ]
        ):
            return False
    return True
