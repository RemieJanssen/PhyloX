def is_endpoint_of_w_fence(self, node):
    if not self.is_reticulation(node):
        return False
    previous_node = node
    current_node = self.child(node)
    currently_at_fence_top = False
    while True:
        if self.is_leaf(current_node):
            return False
        if self.is_reticulation(current_node):
            if currently_at_fence_top:
                return True
            next_node = self.parent(current_node, exclude=[previous_node])
        if self.is_tree_node(current_node):
            if not currently_at_fence_top:
                return False
            next_node = self.child(current_node, exclude=[previous_node])
        previous_node, current_node = current_node, next_node
        currently_at_fence_top = not currently_at_fence_top


def is_tree_based(network):
    if not network.is_binary():
        raise CannotComputeError(
            "tree-basedness cannot be computed for non-binary networks yet."
        )

    if len(network) > 0 and not nx.is_weakly_connected(network):
        return False

    if len(network.roots) > 1:
        return False

    for node in network.nodes:
        if network.is_endpoint_of_w_fence(node):
            return False
    return True
