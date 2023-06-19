import unittest

import pytest

from phylox import DiNetwork
from phylox.exceptions import InvalidMoveDefinitionException, InvalidMoveException
from phylox.rearrangement.move import Move, apply_move, apply_move_sequence
from phylox.rearrangement.movetype import MoveType
from phylox.rearrangement.rearrangementproblem import RearrangementProblem


class TestRearrangementProblem(unittest.TestCase):
    @staticmethod
    def setup_simple_problem():
        nw_1 = DiNetwork(
            nodes=[10],
            edges=[[0, 1], [1, 2], [1, 3], [2, 3], [2, 4], [3, 5]],
            labels=[[4, 1], [5, 2]],
        )
        nw_2 = DiNetwork(
            nodes=[10],
            edges=[[0, 1], [1, 2], [1, 3], [3, 2], [2, 4], [3, 5]],
            labels=[[4, 1], [5, 2]],
        )
        move_type = MoveType.RSPR
        return RearrangementProblem(nw_1, nw_2, move_type)

    def test_simple_problem(self):
        problem = self.setup_simple_problem()
        self.assertEqual(problem.move_type, MoveType.RSPR)

    def test_apply_valid_move(self):
        m = Move(
            origin=(2, 5),
            moving_edge=(1, 3),
            target=(2, 4),
            move_type=MoveType.HEAD,
        )
        problem = self.setup_simple_problem()
        network1_applied = apply_move(problem.network1, m)
        assert True

    def test_check_solution_valid(self):
        problem = self.setup_simple_problem()
        m = Move(
            origin=(2, 5),
            moving_edge=(1, 3),
            target=(2, 4),
            move_type=MoveType.HEAD,
        )
        assert problem.check_solution([m])

    def test_check_solution_valid_partial_isom(self):
        problem = self.setup_simple_problem()
        m = Move(
            origin=(2, 5),
            moving_edge=(1, 3),
            target=(2, 4),
            move_type=MoveType.HEAD,
        )
        assert problem.check_solution([m], isomorphism=[(4, 4), (3, 2)])

    def test_check_solution_valid_invalid_partial_isom_label(self):
        problem = self.setup_simple_problem()
        m = Move(
            origin=(2, 5),
            moving_edge=(1, 3),
            target=(2, 4),
            move_type=MoveType.HEAD,
        )
        assert not problem.check_solution([m], isomorphism=[(4, 5)])

    def test_check_solution_valid_invalid_partial_isom_shape(self):
        problem = self.setup_simple_problem()
        m = Move(
            origin=(2, 5),
            moving_edge=(1, 3),
            target=(2, 4),
            move_type=MoveType.HEAD,
        )
        assert not problem.check_solution([m], isomorphism=[(3, 3)])

    def test_check_solution_wrong_move_type(self):
        problem = self.setup_simple_problem()
        problem.move_type = MoveType.TAIL
        m = Move(
            origin=(2, 5),
            moving_edge=(1, 3),
            target=(2, 4),
            move_type=MoveType.HEAD,
        )
        assert not problem.check_solution([m], isomorphism=[(4, 5), (3, 2)])

    def test_apply_invalid_head_move_cycle(self):
        m = Move(
            origin=(2, 5),
            moving_edge=(1, 3),
            target=(0, 1),
            move_type=MoveType.HEAD,
        )
        problem = self.setup_simple_problem()
        with pytest.raises(
            InvalidMoveException, match="reattachment would create a cycle"
        ):
            network1_applied = apply_move(problem.network1, m)

    def test_apply_invalid_tail_move_cycle(self):
        m = Move(
            origin=(0, 2),
            moving_edge=(1, 3),
            target=(3, 5),
            move_type=MoveType.TAIL,
        )
        problem = self.setup_simple_problem()
        with pytest.raises(
            InvalidMoveException, match="reattachment would create a cycle"
        ):
            network1_applied = apply_move(problem.network1, m)

    def test_apply_invalid_tail_move_not_movable_triangle(self):
        m = Move(
            origin=(1, 3),
            moving_edge=(2, 4),
            target=(0, 1),
            move_type=MoveType.TAIL,
        )
        problem = self.setup_simple_problem()
        with pytest.raises(
            InvalidMoveException, match="removal creates parallel edges"
        ):
            network1_applied = apply_move(problem.network1, m)

    def test_apply_invalid_head_move_not_movable_triangle(self):
        network = DiNetwork(
            edges=[[0, 1], [2, 3], [1, 3], [1, 4], [3, 4], [4, 5]],
        )

        m = Move(
            origin=(1, 4),
            moving_edge=(2, 3),
            target=(4, 5),
            move_type=MoveType.HEAD,
        )
        with pytest.raises(
            InvalidMoveException, match="removal creates parallel edges"
        ):
            apply_move(network, m)

    def test_apply_invalid_bad_origin(self):
        network = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3], [3, 4]],
        )
        m = Move(
            origin=(3, 4),
            moving_edge=(1, 2),
            target=(0, 1),
            move_type=MoveType.TAIL,
        )
        with pytest.raises(
            InvalidMoveException,
            match="origin does not match parent and child or moving_endpoint",
        ):
            apply_move(network, m)

    def test_apply_invalid_tail_move_parallel(self):
        network = DiNetwork(
            edges=[[0, 1], [1, 2], [3, 4], [4, 5], [4, 1]],
        )

        m = Move(
            origin=(3, 5),
            moving_edge=(4, 1),
            target=(0, 1),
            move_type=MoveType.TAIL,
        )
        with pytest.raises(
            InvalidMoveException, match="reattachment creates parallel edges"
        ):
            apply_move(network, m)

    def test_apply_invalid_head_move_parallel(self):
        network = DiNetwork(
            edges=[[0, 1], [1, 2], [3, 4], [4, 5], [4, 1]],
        )

        m = Move(
            origin=(0, 2),
            moving_edge=(4, 1),
            target=(4, 5),
            move_type=MoveType.HEAD,
        )
        with pytest.raises(
            InvalidMoveException, match="reattachment creates parallel edges"
        ):
            apply_move(network, m)

    def test_apply_valid_vertical_plus(self):
        network1 = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3]],
            labels=[[2, "a"], [3, "b"]],
        )
        network2 = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3], [2, 3], [2, 4], [3, 5]],
            labels=[[4, "a"], [5, "b"]],
        )
        m = Move(
            start_edge=(1, 2),
            end_edge=(1, 3),
            network=network1,
            move_type=MoveType.VPLU,
        )
        problem = RearrangementProblem(network1, network2, move_type=MoveType.VPLU)
        problem.check_solution([m])

    def test_apply_invalid_vertical_plus(self):
        network = DiNetwork(
            edges=[[0, 1], [1, 2], [1, 3]],
        )
        m = Move(
            start_edge=(1, 2),
            end_edge=(0, 1),
            network=network,
            move_type=MoveType.VPLU,
        )
        with pytest.raises(
            InvalidMoveException, match="end node is reachable from start node"
        ):
            apply_move(network, m)
