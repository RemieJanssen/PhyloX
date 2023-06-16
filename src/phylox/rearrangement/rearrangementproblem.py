from copy import deepcopy

from phylox.isomorphism import is_isomorphic
from phylox.rearrangement.move import apply_move_sequence


class RearrangementProblem(object):
    def __init__(self, network1, network2, move_type):
        self.network1 = network1
        self.network2 = network2
        self.move_type = move_type

    def check_solution(self, seq_moves, isomorphism=None):
        if not all([move.is_type(self.move_type) for move in seq_moves]):
            return False

        final_network = apply_move_sequence(self.network1, seq_moves)

        return is_isomorphic(
            final_network, self.network2, partial_isomorphism=isomorphism
        )
