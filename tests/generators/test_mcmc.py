import unittest

import pytest

from phylox import DiNetwork
from phylox.classes.dinetwork import is_stack_free
from phylox.generators.mcmc import sample_mcmc_networks
from phylox.isomorphism import is_isomorphic
from phylox.rearrangement.movetype import MoveType

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

    def test_seed(self):
        number_of_samples = 2
        network1 = DiNetwork(
            edges=[(0, 1), (0, 2), (1, 3), (1, 4), (2, 5), (2, 6)],
            labels=[(3, "A"), (4, "B"), (5, "C"), (6, "D")]
        )
        move_type_probabilities = {
            MoveType.TAIL: 0.5,
            MoveType.HEAD: 0.5,
        }
        samples1 = sample_mcmc_networks(
            network1,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=number_of_samples,
            burn_in=10,
            seed=1,
        )
        samples2 = sample_mcmc_networks(
            network1,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=number_of_samples,
            burn_in=10,
            seed=1,
        )
        samples3 = sample_mcmc_networks(
            network1,
            move_type_probabilities=move_type_probabilities,
            number_of_samples=number_of_samples,
            burn_in=10,
            seed=2,
        )
        self.assertEqual(len(samples1), number_of_samples)
        self.assertEqual(len(samples2), number_of_samples)
        self.assertEqual(len(samples3), number_of_samples)
        for sample1, sample2 in zip(samples1, samples2):
            self.assertTrue(is_isomorphic(sample1, sample2))
        for sample1, sample3 in zip(samples1, samples3):
            self.assertFalse(is_isomorphic(sample1, sample3))