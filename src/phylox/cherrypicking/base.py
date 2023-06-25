from copy import deepcopy
from enum import Enum


class CHERRYTYPE(Enum):
    CHERRY = 1
    RETICULATEDCHERRY = 2
    NONE = 0


def find_reducible_pairs_with_second(N, x):
    """
    Finds a list of reducible pairs (cherries and reticulated cherries) in the
    network N with leaf x as second element of the pair.
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
            for sibling_child in N.successors(sibling):
                if N.out_degree(sibling_child) == 0:
                    reducible_pairs.append((sibling_child, x))
    return reducible_pairs


def find_reticulated_cherry_with_first(N, x):
    """
    Finds a list of reticulated cherries in the network N with leaf x as first
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


def reduce_pair(network, x, y):
    network = deepcopy(network)

    cherry_type = check_reducible_pair(network, x, y)
    if cherry_type == CHERRYTYPE.CHERRY:
        px = network.parent(x)
        network.remove_node(x)
        if network.out_degree(px) == 1:
            ppx = network.parent(px)
            network.remove_node(px)
            network.add_edge(ppx, y)
    if cherry_type == CHERRYTYPE.RETICULATEDCHERRY:
        px = network.parent(x)
        py = network.parent(y)
        network.remove_edge(py, px)
        if network.in_degree(px) == 1:
            ppx = network.parent(px)
            network.add_edge(ppx, x)
            network.remove_node(px)
        if network.out_degree(py) == 1:
            ppy = network.parent(py)
            network.add_edge(ppy, y)
            network.remove_node(py)
    return network, cherry_type


def check_reducible_pair(network, x, y):
    if network.has_node(x):
        if network.has_node(y):
            for px in network.predecessors(x):
                for py in network.predecessors(y):
                    if px == py:
                        return CHERRYTYPE.CHERRY
                    if network.out_degree(px) == 1:
                        if px in network.successors(py):
                            return CHERRYTYPE.RETICULATEDCHERRY
    return CHERRYTYPE.NONE
