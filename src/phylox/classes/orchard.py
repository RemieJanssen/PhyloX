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
            pair = network_copy.is_second_in_reducible_pair(leaf)
            print("pair", pair)
            if pair:
                reduced = network_copy.reduce_pair(pair)
                if reduced == "cherry":
                    leaves.remove(pair[0])
                checked_all_leaves = False
                break
        if len(network_copy.edges) == 1:
            return True
        done = checked_all_leaves
    return False