from enum import Enum
from copy import deepcopy

class CHERRYTYPE(Enum):
    CHERRY = 1
    RETICULATEDCHERRY = 2
    NONE = 0

def is_second_in_reducible_pair(network, x):
    for node in network.predecessors(x):
        px = node
    for cpx in network.successors(px):
        if cpx != x:
            if network.out_degree(cpx)==0:
                return (cpx,x)
            if network.out_degree(cpx)==1:
                for ccpx in network.successors(cpx):
                    if network.out_degree(ccpx)==0:
                        return (ccpx,x)
    return False

def reduce_pair(network, x, y):
    network = deepcopy(network)

    cherry_type = check_cherry(network, x, y)
    if cherry_type == CHERRYTYPE.CHERRY:
        px = network.parent(x)
        network.remove_node(x)
        if network.out_degree(px) == 1:
            ppx = network.parent(px)
            network.remove_node(px)
            network.add_edge(ppx,y)
    if cherry_type == CHERRYTYPE.RETICULATEDCHERRY:
        px = network.parent(x)
        py = network.parent(y)
        network.remove_edge(py,px)
        if network.in_degree(px) == 1:
            ppx = network.parent(px)
            network.add_edge(ppx, x)
            network.remove_node(px)
        if network.out_degree(py) == 1:
            ppy = network.parent(py)
            network.add_edge(ppy, y)
            network.remove_node(py)
    return network, cherry_type

def check_cherry(network, x, y):
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
