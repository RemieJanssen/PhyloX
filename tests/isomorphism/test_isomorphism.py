import unittest

from phylox import DiNetwork, LABEL_ATTR
from phylox.isomorphism import is_isomorphic

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
        self.assertFalse(is_isomorphic(network1, network2, partial_isomorphism=[(3, 4), (4, 3)]))

    def test_isomorphism_with_partial_isomorphism_on_leaves_swapped(self):
        network1 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        network2 = DiNetwork(
            edges=[(1, 2), (2, 3), (2, 4), (3, 4), (3, 5), (4, 6)],
            labels=[(5, "a"), (6, "b")],
        )
        self.assertFalse(is_isomorphic(network1, network2, partial_isomorphism=[(6, 5), (5, 6)]))

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
        self.assertFalse(is_isomorphic(network1, network2, ignore_labels=True, partial_isomorphism=[(6, 5), (5, 6)]))
