import unittest
import pytest

from phylox.generators.mcmc import sample_mcmc_networks
from phylox import DiNetwork
from phylox.rearrangement.movetype import MoveType
from phylox.isomorphism import is_isomorphic
from phylox.classes.dinetwork import is_stack_free


class TestMCMCSamples(unittest.TestCase):
    def test_sample_no_moves(self):
        network = DiNetwork(
            edges=[(0, 1), (0, 2), (1, 3), (1, 4)],
        )
        move_type_probabilities = {
            MoveType.TAIL: 0.4,
            MoveType.HEAD: 0.4,
            MoveType.VPLU: 0.1,
            MoveType.VMIN: 0.1,
        }
        samples = sample_mcmc_networks(
            network,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=1,
            burn_in=0,
        )
        self.assertEqual(len(samples), 1)
        self.assertEqual(len(samples[0].edges()), 4)
        self.assertTrue(is_isomorphic(network, samples[0]))

    def test_sample_horizontal(self):
        network = DiNetwork(
            edges=[(0, 1), (0, 2), (1, 3), (1, 4)],
        )
        move_type_probabilities = {
            MoveType.TAIL: 0.5,
            MoveType.HEAD: 0.5,
        }
        samples = sample_mcmc_networks(
            network,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=10,
            burn_in=10,
        )
        self.assertEqual(len(samples), 10)
        for sample in samples:
            self.assertEqual(len(sample.edges()), 4)
            self.assertTrue(is_isomorphic(network, sample))

    def test_sample_horizontal_add_root(self):
        network = DiNetwork(
            edges=[(0, 1), (0, 2), (1, 3), (1, 4)],
        )
        move_type_probabilities = {
            MoveType.TAIL: 0.5,
            MoveType.HEAD: 0.5,
        }
        samples = sample_mcmc_networks(
            network,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=10,
            burn_in=10,
            add_root_if_necessary=True,
        )
        self.assertEqual(len(samples), 10)
        for sample in samples:
            self.assertEqual(len(sample.edges()), 5)

    def test_sample_vertical(self):
        network = DiNetwork(
            edges=[(0, 1), (0, 2), (1, 3), (1, 4)],
        )
        move_type_probabilities = {
            MoveType.VPLU: 0.5,
            MoveType.VMIN: 0.5,
        }
        samples = sample_mcmc_networks(
            network,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=3,
            burn_in=3,
        )
        self.assertEqual(len(samples), 3)

    def test_sample_vertical_ignore_symmetry(self):
        network = DiNetwork(
            edges=[(0, 1), (0, 2), (1, 3), (1, 4)],
        )
        move_type_probabilities = {
            MoveType.VPLU: 0.5,
            MoveType.VMIN: 0.5,
        }
        samples = sample_mcmc_networks(
            network,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=10,
            burn_in=20,
            correct_symmetries=False,
        )
        self.assertEqual(len(samples), 10)

    def test_sample_horizontal_stack_free(self):
        network = DiNetwork(
            edges=[
                (1, 2),
                (2, 3),
                (2, 4),
                (3, 5),
                (3, 6),
                (4, 5),
                (4, 6),
                (5, 7),
                (6, 8),
            ],
        )
        move_type_probabilities = {
            MoveType.TAIL: 0.5,
            MoveType.HEAD: 0.5,
        }
        samples = sample_mcmc_networks(
            network,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=100,
            burn_in=25,
            correct_symmetries=False,
            restriction_map=is_stack_free,
        )
        self.assertEqual(len(samples), 100)
        for sample in samples:
            self.assertTrue(is_stack_free(sample))
