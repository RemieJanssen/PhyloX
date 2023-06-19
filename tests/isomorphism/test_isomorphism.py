import unittest

from phylox import LABEL_ATTR, DiNetwork
from phylox.isomorphism import count_automorphisms, is_isomorphic


class TestIsomorphism(unittest.TestCase):
    def test_isomorphism(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
        )
        self.assertTrue(is_isomorphic(network1, network2))

    def test_isomorphism_with_labels(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertTrue(is_isomorphic(network1, network2))

    def test_isomorphism_with_switched_labels(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(6, "a"), (5, "b")],
        )
        self.assertFalse(is_isomorphic(network1, network2))

    def test_isomorphism_with_partial_isomorphism(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertTrue(is_isomorphic(network1, network2, partial_isomorphism=[(1, 1)]))

    def test_isomorphism_with_partial_isomorphism_on_leaves(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertTrue(is_isomorphic(network1, network2, partial_isomorphism=[(5, 5)]))

    def test_isomorphism_with_partial_isomorphism_on_leaves_swapped(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertFalse(
            is_isomorphic(network1, network2, partial_isomorphism=[(3, 4), (4, 3)])
        )

    def test_isomorphism_with_partial_isomorphism_on_leaves_swapped2(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertFalse(
            is_isomorphic(network1, network2, partial_isomorphism=[(6, 5), (5, 6)])
        )

    def test_isomorphism_ignore_labels(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertTrue(is_isomorphic(network1, network2, ignore_labels=True))

    def test_isomorphism_ignore_labels_partial_isom(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertFalse(
            is_isomorphic(
                network1,
                network2,
                ignore_labels=True,
                partial_isomorphism=[(6, 5), (5, 6)],
            )
        )


class TestAutomorphism(unittest.TestCase):
    def test_automorphism_simple(self):
        network = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4)],
            labels=[(3, "a"), (4, "b")],
        )
        self.assertEqual(count_automorphisms(network), 1)
        self.assertEqual(count_automorphisms(network, ignore_labels=True), 2)

    def test_automorphism_larger(self):
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
            labels=[(7, "a"), (8, "b")],
        )
        self.assertEqual(count_automorphisms(network), 2)
        self.assertEqual(count_automorphisms(network, ignore_labels=True), 4)

    def test_automorphism_larger_strange_labels(self):
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
            labels=[(7, "a"), (8, "a")],
        )
        self.assertEqual(count_automorphisms(network), 4)
        self.assertEqual(count_automorphisms(network, ignore_labels=True), 4)
