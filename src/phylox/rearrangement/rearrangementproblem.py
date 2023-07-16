from copy import deepcopy

from phylox.isomorphism import is_isomorphic
from phylox.rearrangement.move import apply_move_sequence


class RearrangementProblem(object):
    """
    A rearrangement problem is a tuple (N1, N2, M) where N1 and N2 are phylogenetic networks and M is a move type.

    :param network1: a phylogenetic network phylox.DiNetwork.
    :param network2: a phylogenetic network phylox.DiNetwork.
    :param move_type: a move type phylox.rearrangement.move.MoveType.

    :example:
    >>> from phylox import DiNetwork
    >>> from phylox.rearrangement.rearrangementproblem import RearrangementProblem
    >>> from phylox.rearrangement.move import MoveType, Move
    >>> network1 = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3)],
    ...     labels=[(2, "A"), (3, "B")],
    ... )
    >>> network2 = DiNetwork(
    ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
    ...     labels=[(4, "A"), (5, "B")],
    ... )
    >>> problem = RearrangementProblem(network1, network2, MoveType.ALL)
    >>> problem.check_solution([
    ...     Move(
    ...         move_type=MoveType.VPLU,
    ...         start_edge=(1,2),
    ...         end_edge=(1,3),
    ...         network = network1,
    ...     ),
    ... ])
    True
    """

    def __init__(self, network1, network2, move_type):
        self.network1 = network1
        self.network2 = network2
        self.move_type = move_type

    def check_solution(self, seq_moves, isomorphism=None):
        """
        Checks if a sequence of moves solves the rearrangement problem.

        :param seq_moves: a sequence of moves phylox.rearrangement.move.Move.
        :param isomorphism: a partial isomorphism between the networks.
        :return: true if the sequence of moves solves the rearrangement problem, false otherwise.

        :example:
        >>> from phylox import DiNetwork
        >>> from phylox.rearrangement.rearrangementproblem import RearrangementProblem
        >>> from phylox.rearrangement.move import MoveType, Move
        >>> network1 = DiNetwork(
        ...     edges=[(0,1),(1,2),(1,3)],
        ...     labels=[(2, "A"), (3, "B")],
        ... )
        >>> network2 = DiNetwork(
        ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
        ...     labels=[(4, "A"), (5, "B")],
        ... )
        >>> problem = RearrangementProblem(network1, network2, MoveType.ALL)
        >>> problem.check_solution([
        ...     Move(
        ...         move_type=MoveType.VPLU,
        ...         start_edge=(1,2),
        ...         end_edge=(1,3),
        ...         network = network1,
        ...     ),
        ... ])
        True
        """
        if not all([move.is_type(self.move_type) for move in seq_moves]):
            return False

        final_network = apply_move_sequence(self.network1, seq_moves)

        return is_isomorphic(
            final_network, self.network2, partial_isomorphism=isomorphism
        )
