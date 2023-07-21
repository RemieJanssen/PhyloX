# Find the original location of the moving_endpoint,
# That is, the edge from which we remove it.
def from_edge(network, moving_edge, moving_endpoint):
    """
    Finds the original location (from-edge) of the moving_endpoint of an edge if it is moved.

    :param network: a phylogenetic network.
    :param moving_edge: an edge of the network.
    :param moving_endpoint: a node of the network, which must be an endpoint of the edge.
    :return: a pair of nodes (p,c) where p and c are a parent and child of the moving_endpoint such that they are both not part of the moving_edge.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.invertsequence import from_edge
    >>> network = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3)],
    ... )
    >>> from_edge(network, (1,2), 1)
    (0, 3)
    """
    other_parent = network.parent(moving_endpoint, exclude=moving_edge)
    other_child = network.child(moving_endpoint, exclude=moving_edge)
    return (other_parent, other_child)


# def InvertMoveSequence(seq):
#     """
#     Inverts a sequence of moves.

#     :param seq: a list of moves, where each moves is in the format (moving_edge,moving_endpoint,from_edge,to_edge).
#     :return: the inverse list of moves.
#     """
#     newSeq = []
#     for move in reversed(seq):
#         moving_edge, moving_endpoint, from_edge, to_edge = move
#         newSeq.append((moving_edge, moving_endpoint, to_edge, from_edge))
#     return newSeq


# def ReplaceNodeNamesInMoveSequence(seq, isom):
#     """
#     Renames the nodes in a sequence of moves using an isomorphism mapping between two networks.

#     :param seq: a list of moves, implicitly using the nodes of a network.
#     :param isom: a dictionary, containing a bijective mapping from nodes of the networks to another set.
#     :return: a list of moves where the nodes are replaced by their image under the isom mapping.
#     """
#     if type(seq) == int:
#         return isom[seq]
#     return list(map(lambda x: ReplaceNodeNamesInMoveSequence(x, isom), seq))


# def ReplaceNodeNamesByOriginal(network, seq):
#     """
#     Renames the nodes in a sequence of moves by their original names as given in the input. These are stored as a node attribute 'original'.

#     :param network: a phylogenetic network.
#     :param seq: a sequence of moves on teh network.
#     :return: a list of moves on the network using the original node names as used in the input.
#     """
#     if type(seq) == int:
#         return network.node[seq]['original']
#     if seq == 'rho':
#         return "rho"
#     return list(map(lambda x: ReplaceNodeNamesByOriginal(network, x), seq))
