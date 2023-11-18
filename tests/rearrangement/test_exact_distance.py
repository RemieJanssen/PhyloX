import unittest

import pytest

from phylox import DiNetwork
from phylox.rearrangement.movetype import MoveType
from phylox.rearrangement.rearrangementproblem import RearrangementProblem


class TestRearrangementProblemExactDistance(unittest.TestCase):
    @staticmethod
    def setup_simple_problem():
        nw_1 = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3], [2, 3], [2, 4], [3, 5]],
            labels=[[4, 1], [5, 2]],
        )
        nw_2 = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3]],
            labels=[[2, 1], [3, 2]],
        )
        move_type = MoveType.ALL
        return RearrangementProblem(nw_1, nw_2, move_type)

    def test_solve_depth_first(self):
        problem = self.setup_simple_problem()
        solution = problem.solve_depth_first()
        assert len(solution) == 1
        assert solution[0].move_type == MoveType.VMIN
