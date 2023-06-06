def is_stack_free(network):
    for node in network.nodes:
        if network.is_reticulation(node) and any(
            [
                network.is_reticulation(child)
                for child in network.successors(node)
            ]
        ):
            return False
    return True
