"""
A module containing the Green Line heuristic and its random version.
The Green Line heuristic is a heuristic to find a sequence of rSPR or tail moves between two phylogenetic networks.
The algorithms referred to in the comments are from the R Janssen's PhD thesis, "Rearranging Phylogenetic Networks", 2021.
"""

from phylox.rearrangement.heuristics.utils import *
from phylox.rearrangement.movability import check_valid, check_movable
from phylox.rearrangement.move import Move, apply_move, MoveType, apply_move_sequence
from phylox.exceptions import InvalidMoveException, InvalidMoveDefinitionException
from networkx.utils.decorators import py_random_state
from phylox.classes.dinetwork import is_binary, is_leaf_labeled_single_root_network


@py_random_state("seed")
def _GL_Case1_rSPR(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False, seed=None):
    """
    An implementation of Algorithm 2. Finds a sequence of rSPR moves that makes it possible to add the lowest reticulation up to the down-closed isomrophism.

    :param N: a phylogenetic network.
    :param Np: a phylogenetic network.
    :param up: a lowest reticulation node of Np above the isomorphism.
    :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
    :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
    :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
    :param seed: a seed for the random number generator.
    :return: a list of rSPR moves in N, a list of rSPR moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism.
    """
    # use notation as in the paper
    # ' is denoted p
    xp = Np.child(up)
    x = isom_Np_N[xp]
    z = N.parent(x, exclude=isom_N_Np.keys(), randomNodes=randomNodes, seed=seed)

    # Case1a: z is a reticulation
    if N.in_degree(z) == 2:
        return [], [], z, up
    # Case1b: z is not a reticulation
    # Find a retic v in N not in the isom yet
    v = FindRetic(N, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes, seed=seed)
    u = None
    for parent in N.predecessors(v):
        try:
            move = Move(
                move_type=MoveType.HEAD,
                moving_edge=(parent, v),
                target=(z, x),
                network=N,
            )
            check_valid(N, move)
        except (InvalidMoveException, InvalidMoveDefinitionException):
            continue
        if not randomNodes:
            u = parent
            break
        elif u == None or seed.getrandbits(1):
            u = parent
    # If v has an incoming arc (u,v) movable to (z,x), perform this move
    if u != None:
        move = Move(
            move_type=MoveType.HEAD, moving_edge=(u, v), target=(z, x), network=N
        )
        return [move], [], v, up
    # if none of the moves is valid, this must be because
    # v is already a reticulation above x with its movable incoming arc (u,v) with u=z.
    return [], [], v, up


@py_random_state("seed")
def GL_Case1_Tail(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False, seed=None):
    """
    An implementation of Algorithm 6. Finds a sequence of tail moves that makes it possible to add the lowest reticulation up to the down-closed isomrophism.

    :param N: a phylogenetic network.
    :param Np: a phylogenetic network.
    :param up: a lowest reticulation node of Np above the isomorphism.
    :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
    :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
    :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
    :param seed: a seed for the random number generator.
    :return: a list of tail moves in N, a list of tail moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism. Returns false if the networks are not isomorphic with 2 leaves and 1 reticulation.
    """
    # use notation as in the paper
    # ' is denoted p
    xp = Np.child(up)
    x = isom_Np_N[xp]
    z = N.parent(x, exclude=isom_N_Np.keys(), randomNodes=randomNodes, seed=seed)
    # Case1a: z is a reticulation
    if N.in_degree(z) == 2:
        return [], [], z, up
    # Case1b: z is not a reticulation
    # z is a tree node
    if check_movable(N, (z, x), z):
        # Case1bi: (z,x) is movable
        # Find a reticulation u in N not in the isomorphism yet
        # TODO: Can first check if the other parent of x suffices here, should heuristcally be better
        u = FindRetic(N, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes)
        v = N.child(u)
        if v == x:
            return [], [], u, up
        # we may now assume v!=x
        if z == v:
            v = N.child(z, exclude=[x])
            w = N.parent(u, randomNodes=randomNodes, seed=seed)
            move = Move(
                move_type=MoveType.TAIL,
                moving_edge=(z, v),
                target=(w, u),
                network=N,
            )
            return [move], [], u, up
        w = Parent(N, u, exclude=[z], randomNodes=randomNodes)
        move1 = Move(
            move_type=MoveType.TAIL, moving_edge=(z, x), target=(u, v), network=N
        )
        move2 = Move(
            move_type=MoveType.TAIL, moving_edge=(z, v), target=(w, u), origin=(u, x)
        )
        return [(move1, move2)], [], u, up
    # Case1bii: (z,x) is not movable
    c = N.parent(z)
    d = N.child(z, exclude=[x])
    # TODO: b does not have to exist if we have an outdeg-2 root, this could be c!
    b = N.parent(c)
    if N.in_degree(b) != 0:
        # Case1biiA: b is not the root of N
        a = N.parent(b, randomNodes=randomNodes, seed=seed)
        # First do the move ((c,d),c,(a,b)), then Case1bi applies as (z,x) is now movable
        move1 = Move(
            move_type=MoveType.TAIL, moving_edge=(c, d), target=(a, b), network=N
        )
        newN = apply_move(N, move1)
        u = FindRetic(
            newN, excludedSet=isom_N_Np.keys(), randomNodes=randomNodes, seed=seed
        )
        v = newN.child(u)
        if v == x:
            # In this case, u is a reticulation parent of x and u is not in the isom. Hence, we can simply add it to the isom.
            # Note: The tail move we did is not necessary!
            # TODO: First check both parents of x for being a reticulation not in the isomorphism yet
            return [], [], u, up
        # we may now assume v!=x
        if z == v:
            # This only happens if z==v and u==b
            # we move z up above the retic b as well, too
            w = newN.parent(b, randomNodes=randomNodes, seed=seed)
            move2 = Move(
                move_type=MoveType.TAIL, moving_edge=(z, d), target=(w, b), network=newN
            )
            return [move1, move2], [], u, up
        w = newN.parent(u, exclude=[z], randomNodes=randomNodes, seed=seed)
        move2 = Move(
            move_type=MoveType.TAIL, moving_edge=(z, x), target=(u, v), network=newN
        )
        move3 = Move(
            move_type=MoveType.TAIL, moving_edge=(z, v), target=(w, u), origin=(u, x)
        )
        return [move1, move2, move3], [], u, up
    # Case1biiB: b is the root of N
    # Note: d is not in the isomorphism
    e = N.child(d)
    if e == x:
        return [], [], d, up
    if N.out_degree(x) == 2:
        s = N.child(x)
        t = N.child(x, exclude=[s])
        if s == e:
            move = Move(
                move_type=MoveType.TAIL, moving_edge=(x, t), target=(d, e), network=N
            )
            return [move], [], d, up
        if t == e:
            move = Move(
                move_type=MoveType.TAIL, moving_edge=(x, s), target=(d, e), network=N
            )
            return [move], [], d, up
        moves = [
            Move(move_type=MoveType.TAIL, moving_edge=(x, s), target=(d, e), network=N),
            Move(
                move_type=MoveType.TAIL,
                moving_edge=(x, e),
                target=(z, t),
                origin=(d, s),
            ),
            Move(
                move_type=MoveType.TAIL,
                moving_edge=(x, t),
                target=(d, s),
                origin=(z, e),
            ),
        ]
        return moves, [], d, up
    if N.out_degree(e) == 2:
        s = N.child(e)
        t = N.child(e, exclude=[s])
        if s == x:
            move = Move(
                move_type=MoveType.TAIL, moving_edge=(e, t), target=(z, x), network=N
            )
            return [move], [], d, up
        if t == x:
            move = Move(
                move_type=MoveType.TAIL, moving_edge=(e, s), target=(z, x), network=N
            )
            return [move], [], d, up
        moves = [
            Move(move_type=MoveType.TAIL, moving_edge=(e, s), target=(z, x), network=N),
            Move(
                move_type=MoveType.TAIL,
                moving_edge=(e, x),
                target=(d, t),
                origin=(z, s),
            ),
            Move(
                move_type=MoveType.TAIL,
                moving_edge=(e, t),
                target=(z, s),
                origin=(d, x),
            ),
        ]
        return moves, [], d, up
    # neither are tree nodes, so both must be leaves
    # In that case, there is no sequence between the two networks.
    return [], [], None, None


@py_random_state("seed")
def GL_Case3(N, Np, up, isom_N_Np, isom_Np_N, randomNodes=False, seed=None):
    """
    An implementation of Algorithm 3. Finds a sequence of tail moves that makes it possible to add the lowest tree node up to the down-closed isomrophism.

    :param N: a phylogenetic network.
    :param Np: a phylogenetic network.
    :param up: a lowest tree node of Np above the isomorphism.
    :param isom_N_Np: a dictionary, containing a partial (down-closed) isomorphism map from N to Np. The inverse of isom_Np_N.
    :param isom_Np_N: a dictionary, containing a partial (down-closed) isomorphism map from Np to N. The inverse of isom_N_Np.
    :param randomNodes: a boolean value, determining whether the random version of this lemma is used.
    :param seed: a seed for the random number generator.
    :return: a list of tail moves in N, a list of tail moves in Np, a node of N, a node of Np. After performing the lists of moves on the networks, the nodes can be added to the isomorphism.
    """
    # Find the children x' and y' of u'
    xp, yp = list(Np.successors(up))
    # Find the corresponding nodes x and y in N
    x = isom_Np_N[xp]
    y = isom_Np_N[yp]
    # Find the set of common parents of x and y
    parents_x = set(N.predecessors(x))
    parents_y = set(N.predecessors(y))
    common_parents = parents_x & parents_y
    # Case3a: x and y have a common parent not in the isom
    common_parents_not_Y = []
    for parent in common_parents:
        if parent not in isom_N_Np.keys():
            # then parent can be mapped to up
            common_parents_not_Y += [parent]
            if not randomNodes:
                return [], [], parent, up
    if common_parents_not_Y:
        return [], [], seed.choice(common_parents_not_Y), up

    # Case3b: x and y do not have a common parent in the isomorphism
    # For both, find a parent not in the isomorphism yet
    # TODO: preferably make them tree nodes
    z_x = N.parent(x, exclude=isom_N_Np.keys(), randomNodes=randomNodes, seed=seed)
    z_y = N.parent(y, exclude=isom_N_Np.keys(), randomNodes=randomNodes, seed=seed)

    # Case3bi: (z_x,x) is movable
    try:
        move = Move(
            move_type=MoveType.TAIL, moving_edge=(z_x, x), target=(z_y, y), network=N
        )
        check_valid(N, move)
        return [move], [], z_x, up
    except (InvalidMoveException, InvalidMoveDefinitionException):
        pass
    # Case3bii: (z_y,y) is movable
    try:
        move = Move(
            move_type=MoveType.TAIL, moving_edge=(z_y, y), target=(z_x, x), network=N
        )
        check_valid(N, move)
        return [move], [], z_y, up
    except (InvalidMoveException, InvalidMoveDefinitionException):
        pass
    # Case3biii: Neither (z_x,x) nor (z_y,y) is movable

    if N.in_degree(z_x) == 2 or N.in_degree(z_y) == 2:
        return [], [], None, None
    # Both these parents are tree nodes.
    # This happens always in the non-random version, as otherwise there is a lowest reticulation node and we would be in Case 1 or 2.

    # As both nodes are tree nodes and the arcs immovable, both arcs hang of the side of a triangle.
    # Find the top node of the triangle for z_x
    c_x = N.parent(z_x)
    b_x = N.parent(c_x)

    # Find the top node of the triangle for z_y
    c_y = N.parent(z_y)
    b_y = N.parent(c_y)

    if N.in_degree(b_x) == 0:
        # c_x is the child of the root
        # c_x!=c_y, so c_y is not the child of the root
        # swap the roles of x and y
        x, y = y, x
        z_x, z_y = z_y, z_x
        b_x, b_y = b_y, b_x
        c_x, c_y = c_y, c_x
    # c_x is not the child of the root
    # find a parent of b_x, and the bottom node of the triangle d_x
    a_x = N.parent(b_x, randomNodes=randomNodes, seed=seed)
    d_x = N.child(c_x, exclude=[z_x])
    # e_x not in the proof, but needed to define the origin of the second move
    e_x = N.parent(c_x)
    moves = [
        Move(
            move_type=MoveType.TAIL,
            moving_edge=(c_x, d_x),
            target=(a_x, b_x),
            network=N,
        ),
        Move(
            move_type=MoveType.TAIL,
            moving_edge=(z_x, x),
            target=(z_y, y),
            origin=(e_x, d_x),
        ),
    ]
    return moves, [], z_x, up


class HeuristicDistanceMixin:
    """
    A class containing the Green Line heuristic and its random version.
    Meant to be inherited by the RearrangementProblem class.
    """

    def check_green_line_requirements(self):
        if not self.move_type in [MoveType.TAIL, MoveType.RSPR, MoveType.ALL]:
            raise Exception("Move type not supported by Green Line heuristic")
        if not is_binary(self.network1) or not is_binary(self.network2):
            raise Exception("Green Line heuristic only works for binary networks")
        if not is_leaf_labeled_single_root_network(
            self.network1
        ) or not is_leaf_labeled_single_root_network(self.network2):
            raise Exception(
                "Green Line heuristic only works for leaf-labeled networks with a single root"
            )
        if not len(self.network1.leaves) == len(self.network1.labels) or not len(
            self.network2.leaves
        ) == len(self.network2.labels):
            raise Exception(
                "Green Line heuristic only works for networks with unique labels"
            )
        if self.network1.labels == self.network2.labels:
            raise Exception(
                "Green Line heuristic only works for networks with the same set of labels"
            )

    def heuristic_green_line(self):
        """
        An implementation of Algorithm 4 and its tail move counterpart. Finds a sequence of tail/rSPR moves from network1 to network2 by building a down-closed isomorphism.
        Assumes the networks have the same leaf set, the same number of reticulations, are both binary, and all labels are unique.

        :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
        """
        self.check_green_line_requirements()

        head_moves = self.move_type in [MoveType.RSPR, MoveType.ALL]

        # Find the root and labels of the networks
        root1 = list(self.network1.roots)[0]
        root2 = list(self.network2.roots)[0]

        # initialize isomorphism
        isom_1_2 = dict()
        isom_2_1 = dict()
        isom_size = 0
        for label, [node1] in self.network1.labels.items():
            node2 = self.network2.labels[label][0]
            isom_1_2[node1] = node2
            isom_2_1[node2] = node1
            isom_size += 1

        # Keep track of the size of the isomorphism and the size it is at the end of the green line algorithm
        goal_size = len(self.network1) - 1

        # init lists of sequence of moves
        # list of (moving_edge,moving_endpoint,from_edge,to_edge)
        seq_from_1 = []
        seq_from_2 = []
        # TODO keep track of lowest nodes? (Currently, the code does not do this, but it could speed up the code)

        network1 = self.network1
        network2 = self.network2

        # Do the green line algorithm
        while isom_size < goal_size:
            # Find lowest nodes above the isom in the networks:
            (
                lowest_tree_node_network1,
                lowest_retic_network1,
            ) = LowestReticAndTreeNodeAbove(network1, isom_1_2.keys())
            (
                lowest_tree_node_network2,
                lowest_retic_network2,
            ) = LowestReticAndTreeNodeAbove(network2, isom_2_1.keys())

            ######################################
            # Case1: a lowest retic in network1
            if lowest_retic_network1 != None:
                # use notation as in the paper network1 = N', network2 = N, where ' is denoted p
                up = lowest_retic_network1
                if head_moves:
                    (
                        moves_network_2,
                        moves_network_1,
                        added_node_network_2,
                        added_node_network_1,
                    ) = _GL_Case1_rSPR(network2, network1, up, isom_2_1, isom_1_2)
                else:
                    (
                        moves_network_2,
                        moves_network_1,
                        added_node_network_2,
                        added_node_network_1,
                    ) = GL_Case1_Tail(network2, network1, up, isom_2_1, isom_1_2)
                    if added_node_network_1 == None:
                        return False
            ######################################
            # Case2: a lowest retic in network2
            elif lowest_retic_network2 != None:
                # use notation as in the paper network2 = N', network1 = N, where ' is denoted p
                up = lowest_retic_network2
                if head_moves:
                    (
                        moves_network_1,
                        moves_network_2,
                        added_node_network_1,
                        added_node_network_2,
                    ) = _GL_Case1_rSPR(network1, network2, up, isom_1_2, isom_2_1)
                else:
                    (
                        moves_network_1,
                        moves_network_2,
                        added_node_network_1,
                        added_node_network_2,
                    ) = GL_Case1_Tail(network1, network2, up, isom_1_2, isom_2_1)
                    if added_node_network_1 == None:
                        return False

                        ######################################
            # Case3: a lowest tree node in network1
            else:
                # use notation as in the paper network1 = N, network2 = N'
                up = lowest_tree_node_network2
                (
                    moves_network_1,
                    moves_network_2,
                    added_node_network_1,
                    added_node_network_2,
                ) = GL_Case3(network1, network2, up, isom_1_2, isom_2_1)
            # Now perform the moves and update the isomorphism
            isom_1_2[added_node_network_1] = added_node_network_2
            isom_2_1[added_node_network_2] = added_node_network_1
            seq_from_1 += moves_network_1
            seq_from_2 += moves_network_2
            network1 = apply_move_sequence(network1, moves_network_1)
            network2 = apply_move_sequence(network2, moves_network_2)
            isom_size += 1

        # Add the root to the isomorphism, if it was there
        isom_1_2[root1] = root2
        isom_2_1[root2] = root1
        # invert seq_from_2, rename to node names of network1, and append to seq_from_1
        return seq_from_1 + [
            move.invert().rename_nodes(isom_2_1) for move in reversed(seq_from_2)
        ]

    @py_random_state("seed")
    def heuristic_green_line_random(self, seed=None):
        """
        An implementation of Algorithm 5 and its tail move counterpart. Finds a sequence of tail/rSPR moves from network1 to network2 by randomly building a down-closed isomorphism.

        :return: A list of tail/rSPR moves from network1 to network2. Returns False if such a sequence does not exist.
        """
        self.check_green_line_requirements()
        head_moves = self.move_type in [MoveType.RSPR, MoveType.ALL]

        # Find the root and labels of the networks
        root1 = list(self.network1.roots)[0]
        root2 = list(self.network2.roots)[0]

        # initialize isomorphism
        isom_1_2 = dict()
        isom_2_1 = dict()
        isom_size = 0
        for label, [node1] in self.network1.labels.items():
            node2 = self.network2.labels[label][0]
            isom_1_2[node1] = node2
            isom_2_1[node2] = node1
            isom_size += 1

        # Keep track of the size of the isomorphism and the size it is at the end of the green line algorithm
        goal_size = len(self.network1) - 1

        # init lists of sequence of moves
        # list of (moving_edge,moving_endpoint,from_edge,to_edge)
        seq_from_1 = []
        seq_from_2 = []
        # TODO keep track of lowest nodes? (Currently, the code does not do this, but it could speed up the code)

        network1 = self.network1
        network2 = self.network2

        # Do the green line algorithm
        while isom_size < goal_size:
            # Find all lowest nodes above the isom in the networks:
            (
                lowest_tree_node_network1,
                lowest_retic_network1,
            ) = LowestReticAndTreeNodeAbove(network1, isom_1_2.keys(), allnodes=True)
            (
                lowest_tree_node_network2,
                lowest_retic_network2,
            ) = LowestReticAndTreeNodeAbove(network2, isom_2_1.keys(), allnodes=True)

            # Construct a list of all lowest nodes in a tuple with the corresponding network (in random order)
            # I.e. If u is a lowest node of network one, it will appear in the list as (u,1)
            lowest_nodes_network1 = [
                (u, 1) for u in lowest_tree_node_network1 + lowest_retic_network1
            ]
            lowest_nodes_network2 = [
                (u, 2) for u in lowest_tree_node_network2 + lowest_retic_network2
            ]
            candidate_lowest_nodes = lowest_nodes_network1 + lowest_nodes_network2
            seed.shuffle(candidate_lowest_nodes)

            # As some cases do not give an addition to the isom, we continue trying lowest nodes until we find one that does.
            for lowest_node, network_number in candidate_lowest_nodes:
                ######################################
                # Case1: a lowest retic in network1
                if network_number == 1 and network1.in_degree(lowest_node) == 2:
                    # use notation as in the paper network1 = N', network2 = N, where ' is denoted p
                    up = lowest_node
                    if head_moves:
                        (
                            moves_network_2,
                            moves_network_1,
                            added_node_network_2,
                            added_node_network_1,
                        ) = _GL_Case1_rSPR(
                            network2,
                            network1,
                            up,
                            isom_2_1,
                            isom_1_2,
                            randomNodes=True,
                            seed=seed,
                        )
                    else:
                        (
                            moves_network_2,
                            moves_network_1,
                            added_node_network_2,
                            added_node_network_1,
                        ) = GL_Case1_Tail(
                            network2,
                            network1,
                            up,
                            isom_2_1,
                            isom_1_2,
                            randomNodes=True,
                            seed=seed,
                        )
                        if added_node_network_1 == None:
                            # The networks are non-isom networks with 2 leaves and 1 reticulation
                            return False
                    # This case always gives a node to add to the isom
                    break

                ######################################
                # Case2: a lowest retic in network2
                elif network_number == 2 and network2.in_degree(lowest_node) == 2:
                    # use notation as in the paper network2 = N', network1 = N, where ' is denoted p
                    up = lowest_node
                    if head_moves:
                        (
                            moves_network_1,
                            moves_network_2,
                            added_node_network_1,
                            added_node_network_2,
                        ) = _GL_Case1_rSPR(
                            network1,
                            network2,
                            up,
                            isom_1_2,
                            isom_2_1,
                            randomNodes=True,
                            seed=seed,
                        )
                    else:
                        (
                            moves_network_1,
                            moves_network_2,
                            added_node_network_1,
                            added_node_network_2,
                        ) = GL_Case1_Tail(
                            network1,
                            network2,
                            up,
                            isom_1_2,
                            isom_2_1,
                            randomNodes=True,
                            seed=seed,
                        )
                        if added_node_network_1 == None:
                            # The networks are non-isom networks with 2 leaves and 1 reticulation
                            return False
                            # This case always gives a node to add to the isom
                    break

                ######################################
                # Case3: a lowest tree node in network1
                elif network_number == 2 and network2.out_degree(lowest_node) == 2:
                    # use notation as in the paper network1 = N, network2 = N'
                    up = lowest_node
                    (
                        moves_network_1,
                        moves_network_2,
                        added_node_network_1,
                        added_node_network_2,
                    ) = GL_Case3(
                        network1,
                        network2,
                        up,
                        isom_1_2,
                        isom_2_1,
                        randomNodes=True,
                        seed=seed,
                    )
                    # If we can add a node to the isom, added_node_network_2 has a value
                    if added_node_network_2:
                        break

                ######################################
                # Case3': a lowest tree node in network2
                else:
                    # use notation as in the paper network1 = N, network2 = N'
                    up = lowest_node
                    (
                        moves_network_2,
                        moves_network_1,
                        added_node_network_2,
                        added_node_network_1,
                    ) = GL_Case3(
                        network2,
                        network1,
                        up,
                        isom_2_1,
                        isom_1_2,
                        randomNodes=True,
                        seed=seed,
                    )
                    # If we can add a node to the isom, added_node_network_2 has a value
                    if added_node_network_2:
                        break

            # Now perform the moves and update the isomorphism
            isom_1_2[added_node_network_1] = added_node_network_2
            isom_2_1[added_node_network_2] = added_node_network_1
            seq_from_1 += moves_network_1
            seq_from_2 += moves_network_2
            network1 = apply_move_sequence(network1, moves_network_1)
            network2 = apply_move_sequence(network2, moves_network_2)
            isom_size += 1

        # Add the root to the isomorphism, if it was there
        isom_1_2[root1] = root2
        isom_2_1[root2] = root1

        # invert seq_from_2, rename to node names of network1, and append to seq_from_1
        return seq_from_1 + [
            move.invert().rename_nodes(isom_2_1) for move in reversed(seq_from_2)
        ]
