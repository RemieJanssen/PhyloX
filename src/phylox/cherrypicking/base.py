from enum import Enum

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
    cherry_type = check_cherry(network, x, y)
    if cherry_type == CHERRYTYPE.CHERRY:
        for px in network.predecessors(x):
            network.remove_node(x)
            for ppx in network.predecessors(px):
                network.remove_node(px)
                network.add_edge(ppx,y)
            return True
    if cherry_type == CHERRYTYPE.RETICULATEDCHERRY:
        for px in network.predecessors(x):
            for py in network.predecessors(y):
                network.remove_edge(py,px)
                if network.in_degree(px) == 1:
                    for ppx in network.predecessors(px):
                        network.add_edge(ppx, x)
                        network.remove_node(px)
                for ppy in network.predecessors(py):
                    network.add_edge(ppy, y)
                    network.remove_node(py)
                return True
    return False    

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
