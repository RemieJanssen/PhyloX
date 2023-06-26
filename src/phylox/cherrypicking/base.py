from copy import deepcopy
from enum import Enum
from phylox import DiNetwork


class CHERRYTYPE(Enum):
    CHERRY = 1
    RETICULATEDCHERRY = 2
    NONE = 0


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


def reduce_pair(network, x, y, inplace=False):
    if not inplace:
        network = deepcopy(network)

    cherry_type = check_cherry(network, x, y)
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
    if inplace:
        # TODO empty cache for network properties?
        pass
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



def add_pair(network, x, y, height = [1,1], inplace=False):
    """
    Adds a pair to the network, using the construction from a cherry-picking sequence
    :param x: first element of the pair
    :param y: second element of the pair
    :param height: height of the pair
    :return: false if y is not yet in the network and the network is not empty
    """
    if not inplace:
        network = deepcopy(network)

    #if the network is empty, create a cherry (x,y)
    if len(network.leaves)==0:
        network.add_edge(0,1, length=0)
        network.add_edge(1,2, length=height[0])
        network.add_edge(1,3, length=height[1])
        network.labels[x]=2
        network.labels[y]=3
        network.leaf_nodes[2]=x
        network.leaf_nodes[3]=y
        network.no_nodes=4
        return True
    #if y is not in the network return false, as there is no way to add the pair and get a phylogenetic network
    if y not in network.leaves:
        return False
    #add the pair to the existing network
    node_y=network.labels[y]
    parent_node_y = -1
    for p in network.nw.predecessors(node_y):
        parent_node_y=p

    #first add all edges around y
    length_incoming_y = network.nw[parent_node_y][node_y]['length']
    no_of_trees_incoming_y = network.nw[parent_node_y][node_y]['no_of_trees']
    height_goal_x = height[0]
    if height[1]<length_incoming_y:
        height_pair_y_real = height[1]
    else:
        height_pair_y_real = length_incoming_y
        height_goal_x+=height[1]-height_pair_y_real


    network.nw.add_edge(node_y,network.no_nodes,no_of_trees=no_of_trees_incoming_y+len(red_trees-current_trees), length=height_pair_y_real)
    network.nw[parent_node_y][node_y]['length'] = length_incoming_y - height_pair_y_real
    network.leaf_nodes.pop(network.labels[y],False)
    network.labels[y]=network.no_nodes
    network.leaf_nodes[network.no_nodes]=y

    #Now also add edges around x
    #x is not yet in the network, so make a cherry (x,y)
    if x not in network.leaves:
        network.nw.add_edge(node_y,network.no_nodes+1,no_of_trees=len(red_trees), length=height_goal_x)
        network.leaves.add(x)
        network.labels[x]=network.no_nodes+1
        network.leaf_nodes[network.no_nodes+1]=x
        network.no_nodes+=2
    #x is already in the network, so create a reticulate cherry (x,y)
    else:
        node_x=network.labels[x]
        for parent in network.nw.predecessors(node_x):
            px = parent
        length_incoming_x = network.nw[px][node_x]['length']
        no_of_trees_incoming_x = network.nw[px][node_x]['no_of_trees']
        #if x is below a reticulation, and the height of the new pair is above the height of this reticulation, add the new hybrid arc to the existing reticulation
        if network.nw.in_degree(px)>1 and length_incoming_x<=height_goal_x:
            network.nw.add_edge(node_y,px,no_of_trees=len(red_trees), length=height_goal_x-length_incoming_x)
            network.nw[px][node_x]['no_of_trees']+=len(red_trees)
            network.no_nodes+=1
        #create a new reticulation vertex above x to attach the hybrid arc to
        else:
            height_pair_x = min(height_goal_x,length_incoming_x)
            network.nw.add_edge(node_y,node_x,no_of_trees=len(red_trees), length=height_goal_x-height_pair_x)
            network.nw.add_edge(node_x,network.no_nodes+1,no_of_trees = no_of_trees_incoming_x+len(red_trees), length = height_pair_x)
            network.nw[px][node_x]['length'] = length_incoming_x - height_pair_x
            network.leaf_nodes.pop(network.labels[x],False)
            network.labels[x]=network.no_nodes+1
            network.leaf_nodes[network.no_nodes+1]=x
            network.no_nodes+=2
    return True


class CherryPickingMixin:
    @classmethod
    def from_cherry_picking_sequence(cls, sequence, heights=None):
        network = DiNetwork()
        heights = heights or [[1,1]]*total_len
        for pair, height in zip(sequence, heights):
            add_pair(network, *pair, height=height, inplace=True)
        return network

