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

    @staticmethod
    def setup_large_problem():
        nw_1 = DiNetwork(
            edges=[
                [0, 1],
                [1, 2],
                [1, 3],
                [2, 3],
                [2, 4],
                [3, 5],
                [5, 6],
                [5, 7],
                [6, 8],
                [7, 8],
                [8, 9],
                [6, 10],
                [7, 11],
            ],
            labels=[[4, 1], [9, 3], [10, 4], [11, 2]],
        )
        nw_2 = DiNetwork(
            edges=[
                [0, 7],
                [1, 2],
                [1, 3],
                [2, 3],
                [2, 4],
                [3, 5],
                [5, 6],
                [5, 8],
                [6, 8],
                [8, 9],
                [6, 10],
                [7, 11],
                [7, 1],
            ],
            labels=[[11, 1], [9, 2], [4, 4], [10, 3]],
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

    def test_solve_green_line_large(self):
        problem = self.setup_large_problem()
        solution = problem.heuristic_green_line()
        assert problem.check_solution(solution)

    def test_solve_green_line_random_large(self):
        problem = self.setup_large_problem()
        solution1 = problem.heuristic_green_line_random(seed=1)
        solution2 = problem.heuristic_green_line_random(seed=1)
        assert problem.check_solution(solution1)
        assert [move.__dict__ for move in solution1] == [
            move.__dict__ for move in solution2
        ]
