from copy import deepcopy
from phylox.cherrypicking import is_second_in_reducible_pair, reduce_pair

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

def is_orchard(network):
    if len(network) == 0:
        return True
    leaves = network.leaves
    root = list(network.roots)[0]

    # make a copy and fix a root edge
    network_copy = deepcopy(network)
    if network_copy.out_degree(root) > 1:
        new_node = -1
        while new_node in network_copy.nodes:
            new_node -= 1
        network_copy.add_edge(new_node, root)

    # try to reduce the network copy
    done = False
    while not done:
        checked_all_leaves = True
        for leaf in leaves:
            print("leaf", leaf)
            pair = is_second_in_reducible_pair(network_copy, leaf)
            print("pair", pair)
            if pair:
                reduced = reduce_pair(network_copy, *pair)
                if reduced == "cherry":
                    leaves.remove(pair[0])
                checked_all_leaves = False
                break
        if len(network_copy.edges) == 1:
            return True
        done = checked_all_leaves
    return False

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
