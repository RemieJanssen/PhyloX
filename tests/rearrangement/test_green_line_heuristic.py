import unittest

import pytest

from phylox import DiNetwork
from phylox.rearrangement.movetype import MoveType
from phylox.rearrangement.rearrangementproblem import RearrangementProblem


class TestRearrangementProblemGreenLine(unittest.TestCase):
    @staticmethod
    def setup_simple_problem():
        nw_1 = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3], [2, 3], [2, 4], [3, 5]],
            labels=[[4, 1], [5, 2]],
        )
        nw_2 = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3], [2, 3], [2, 4], [3, 5]],
            labels=[[5, 1], [4, 2]],
        )
        move_type = MoveType.ALL
        return RearrangementProblem(nw_1, nw_2, move_type)

    def test_solve_green_line(self):
        problem = self.setup_simple_problem()
        solution = problem.heuristic_green_line()
        assert problem.check_solution(solution)


    def test_solve_green_line_random(self):
        problem = self.setup_simple_problem()
        solution = problem.heuristic_green_line_random()
        assert problem.check_solution(solution)
