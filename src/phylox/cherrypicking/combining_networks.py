import random
import time
from copy import deepcopy

import networkx as nx
import numpy as np

from phylox import DiNetwork
from phylox.cherrypicking import (
    CHERRYTYPE,
    add_roots_to_sequence,
    check_reducible_pair,
    cherry_height,
    find_all_reducible_pairs,
    find_reducible_pairs_with_first,
    find_reducible_pairs_with_second,
    find_reticulated_cherry_with_first,
    get_indices_of_reducing_pairs,
    reduce_pair,
)
from phylox.constants import LABEL_ATTR, LENGTH_ATTR

# prefix for harmonized node names
HARMONIZE_NODES_BY_LABEL_PREFIX = "hnbl__"


class HybridizationProblem:
    """
    A class to represent a hybridization problem.
    I.e. a set of phylogenetic networks that need to be combined into a single phylogenetic network.

    :param list_of_networks: a list of phylogenetic networks, each given as a phylox.DiNetwork.
    :param newick_strings: if True, the input trees are given as newick strings, otherwise as phylox.DiNetworks.
    """

    def __init__(self, list_of_networks=None, newick_strings=True):
        # The dictionary of trees
        self.trees = dict()
        # the set of leaf labels of the trees
        self.leaves = set()

        # the current best sequence we have found for this set of trees
        self.best_seq = None
        # the list of reduced trees for each of the pairs in the best sequence
        self.best_red_trees = None

        # the best sequence for the algorithm using lengths as input as well
        self.best_seq_with_lengths = None
        # the sets of reduced trees for each pair in this sequence
        self.best_seq_with_lengths_red_trees = None
        # the height of each pair in this sequence
        self.best_seq_with_lengths_heights = None

        # true if distances are used
        self.distances = True

        # read the input trees in 'newick_strings'
        list_of_networks = list_of_networks or []
        for n in list_of_networks:
            if newick_strings:
                network = DiNetwork(newick=n)
            else:
                network = n
            self.trees[len(self.trees)] = network
            self.distances = self.distances and all(
                [LENGTH_ATTR in edge[2] for edge in network.edges(data=True)]
            )

        # check that the labels are unique in each tree
        # and that all leaves have a label
        # also set the leaves of the problem
        for i, tree in self.trees.items():
            rename_dict = dict()
            leaf_labels = set()
            for l in tree.leaves:
                leaf_label = tree.nodes[l][LABEL_ATTR]
                if leaf_label in leaf_labels:
                    raise ValueError(
                        "The label {} is not unique in tree {}".format(leaf_label, i)
                    )
                rename_dict[l] = HARMONIZE_NODES_BY_LABEL_PREFIX + leaf_label
                leaf_labels.add(leaf_label)
            self.leaves.update(
                [HARMONIZE_NODES_BY_LABEL_PREFIX + l for l in leaf_labels]
            )
            nx.relabel_nodes(tree, rename_dict, copy=False)
            tree._clear_cached()

    # Find new cherry-picking sequences for the trees and update the best found
    def CPSBound(
        self, repeats=1, progress=False, track=False, lengths=False, time_limit=None
    ):
        """
        Finds a cherry-picking sequence for the input networks, and updates the best sequence found so far.

        The method simply keeps picking (random) cherries until all networks are reduced to a single leaf.
        If lengths is True, the method picks the lowest cherry in each step (CPHeuristicLenths).
        If track is True, the method keeps track of the reducible pairs in each step, and reduces the most common pair (CPHeuristicStorePairs).
        Otherwise, the method simply picks a random cherry in each step (CPHeuristic).


        :param repeats: the number of times to run the algorithm.
        :param progress: if True, print progress information.
        :param track: if True, keep track of reducible pairs and reduce the most common pair.
        :param lengths: if True, pick the lowest cherry in each step.
        :param time_limit: if given, the algorithm stops after this many seconds.
        :return: the best sequence found so far.

        :example:
        >>> from phylox import DiNetwork
        >>> from phylox.cherrypicking.combining_networks import HybridizationProblem
        >>> network1 = DiNetwork(
        ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
        ...     labels=[(4, "A"), (5, "B")],
        ... )
        >>> network2 = DiNetwork(
        ...     edges=[(0,1),(1,2),(1,3),(2,3),(2,4),(3,5)],
        ...     labels=[(4, "A"), (5, "B")],
        ... )
        >>> problem = HybridizationProblem(
        ...     list_of_networks=[network1, network2],
        ...     newick_strings=False,
        ... )
        >>> sequence = problem.CPSBound()
        >>> sequence == [('B', 'A'), ('B', 'A')] or sequence == [('B', 'A'), ('A', 'B')]
        True
        """
        # Set the specific heuristic that we use, based on the user input and whether the trees have lengths
        Heuristic = self.CPHeuristic
        if track and not lengths:
            if progress:
                print("Tracking reducible cherries")
            Heuristic = self.CPHeuristicStorePairs
        if lengths:
            if not self.distances:
                raise ValueError("not all trees have branch lengths!")
            if progress:
                print("Picking the lowest cherry")
            Heuristic = self.CPHeuristicLengths
            heights_best = []
        # Initialize the recorded best sequences and corresponding data
        best = None
        red_trees_best = []
        starting_time = time.time()
        # Try as many times as required by the integer 'repeats'
        for i in range(repeats):
            if lengths:
                new, reduced_trees, seq_heights = Heuristic(progress=progress)
                if progress:
                    print("found sequence of length: " + str(len(new)))
            else:
                new, reduced_trees = Heuristic(progress=progress)
                if progress:
                    print("found sequence of length: " + str(len(new)))
                    print(new)
                    print("improving sequence")
                new, reduced_trees = self.Improve_Sequence(
                    new, reduced_trees, progress=progress
                )
                if progress:
                    print("new length = " + str(len(new)))
                    print(new)
            if progress:
                print("adding roots")
            new, reduced_trees = add_roots_to_sequence(new, reduced_trees)
            if lengths:
                for i in range(len(new) - len(seq_heights)):
                    seq_heights += [seq_heights[-1]]
            if progress:
                print("final length = " + str(len(new)))
            if best == None or len(new) < len(best):
                best = new
                red_trees_best = reduced_trees
                if lengths:
                    heights_best = seq_heights
            if progress:
                print("best sequence has length " + str(len(best)))
            if time_limit and time.time() - starting_time > time_limit:
                break
        new_seq = [
            (
                x[len(HARMONIZE_NODES_BY_LABEL_PREFIX) :],
                y[len(HARMONIZE_NODES_BY_LABEL_PREFIX) :],
            )
            for (x, y) in best
        ]
        if lengths:
            if not self.best_seq_with_lengths or len(new_seq) < len(
                self.best_seq_with_lengths
            ):
                self.best_seq_with_lengths = new_seq
                self.best_seq_with_lengths_red_trees = red_trees_best
                self.best_seq_with_lengths_heights = heights_best
            return self.best_seq_with_lengths
        else:
            if not self.best_seq or len(new_seq) < len(self.best_seq):
                self.best_seq = new_seq
                self.best_red_trees = red_trees_best
            return self.best_seq

    # Version of the code that uses minimal memory: recompute reducible pairs when necessary.
    def CPHeuristic(self, progress=False):
        if progress:
            print("Copying all inputs to reduce on")
        # Works in a copy of the input trees, copy_of_inputs, because trees have to be reduced somewhere.
        copy_of_inputs = deepcopy(self)
        if progress:
            print("Done, starting reduction of trees")
        CPS = []
        reduced_trees = []
        candidate_leaves = deepcopy(self.leaves)
        i = 1
        while copy_of_inputs.trees and i < 10:
            if progress:
                print("Sequence has length: " + str(len(CPS)))
                print(str(len(copy_of_inputs.trees)) + " trees left.\n")
                print("Reducing trivial pairs")
                # First reduce trivial cherries
            new_seq, new_red_trees = copy_of_inputs.Reduce_Trivial_Pairs(
                candidate_leaves
            )
            if progress:
                print("done")
            CPS += new_seq
            reduced_trees += new_red_trees
            if len(copy_of_inputs.trees) == 0:
                break
            # Now reduce a random cherry from a random tree
            random_index, random_tree = random.choice(
                list(copy_of_inputs.trees.items())
            )
            list_of_cherries = find_all_reducible_pairs(random_tree)
            random_cherry = random.choice(list(list_of_cherries))
            CPS += [random_cherry]

            reduced_by_random_cherry = copy_of_inputs.Reduce_Pair_In_All(random_cherry)
            reduced_trees += [reduced_by_random_cherry]
            candidate_leaves = set(random_cherry)
            i += 1
        return CPS, reduced_trees

    # Version of the code that uses more memory: stores all reducible pairs.
    # Runs when user toggles -t or --track
    def CPHeuristicStorePairs(self, progress=False):
        if progress:
            print("Copying all inputs to reduce on")
        # Works in a copy of the input trees, copy_of_inputs, because trees have to be reduced somewhere.
        copy_of_inputs = deepcopy(self)
        if progress:
            print("Done")
        CPS = []
        reduced_trees = []
        candidate_leaves = deepcopy(self.leaves)
        # Make dict of reducible pairs
        if progress:
            print("finding all reducible pairs")
        reducible_pairs = self.Find_All_Pairs()
        if progress:
            print("found all reducible pairs")
        while copy_of_inputs.trees:
            if progress:
                print("Sequence has length: " + str(len(CPS)))
                print(str(len(copy_of_inputs.trees)) + " trees left.\n")
                print("Reducing trivial pairs")
            # First reduce trivial cherries
            (
                new_seq,
                new_red_trees,
                reducible_pairs,
            ) = copy_of_inputs.Reduce_Trivial_Pairs_Store_Pairs(
                candidate_leaves, reducible_pairs
            )
            if progress:
                print("done")
            CPS += new_seq
            reduced_trees += new_red_trees
            if len(copy_of_inputs.trees) == 0:
                break

            # Now reduce a random cherry from a random tree
            # EITHER: (Get random tree, then random pair from the tree), just like in CPHeuristic
            random_index, random_tree = random.choice(
                list(copy_of_inputs.trees.items())
            )
            list_of_cherries = find_all_reducible_pairs(random_tree)
            random_cherry = random.choice(list(list_of_cherries))

            # OR: (Get a random reducible pair from all pairs)
            # Note that this would result in a different algorithm than CPHeuristic, so we use the previous option
            #            random_cherry = random.choice(list(reducible_pairs.keys()))

            CPS += [random_cherry]
            # reduce all trees with this pair, this is where the list of reducible_pairs is used
            # using the list makes it faster to find all trees that need to be reduced.
            reduced_by_random_cherry = copy_of_inputs.Reduce_Pair_In_All(
                random_cherry, reducible_pairs=reducible_pairs
            )
            reducible_pairs = copy_of_inputs.Update_Reducible_Pairs(
                reducible_pairs, reduced_by_random_cherry
            )
            reduced_trees += [reduced_by_random_cherry]
            candidate_leaves = set(random_cherry)
        return CPS, reduced_trees

    # Version of the code that always picks the lowest available pair
    # Runs when user toggles -l or --lengths and all edges in the input trees have lengths.
    def CPHeuristicLengths(self, progress=False):
        if progress:
            print("Copying all inputs to reduce on")
        # Works in a copy of the input trees, copy_of_inputs, because trees have to be reduced somewhere.
        copy_of_inputs = deepcopy(self)
        if progress:
            print("Done")
        CPS = []
        reduced_trees = []
        heights_seq = []

        candidate_leaves = deepcopy(self.leaves)
        # Make dict of reducible pairs
        if progress:
            print("finding all reducible pairs")
        reducible_pairs = self.Find_All_Pairs()
        current_heights = (
            dict()
        )  # for each reducible pair: [0] gives height, [1] the number of trees it was computed in.

        if progress:
            print("found all reducible pairs")
        while copy_of_inputs.trees:
            if progress:
                print("Sequence has length: " + str(len(CPS)))
                print(str(len(copy_of_inputs.trees)) + " trees left.\n")
                # First reduce trivial cherries
                print("Reducing trivial pairs")
            (
                new_seq,
                new_red_trees,
                reducible_pairs,
                new_heights_seq,
            ) = copy_of_inputs.Reduce_Trivial_Pairs_Lengths(
                candidate_leaves, reducible_pairs
            )
            if progress:
                print("done")
            CPS += new_seq
            reduced_trees += new_red_trees
            heights_seq += new_heights_seq
            if len(copy_of_inputs.trees) == 0:
                break

            # Now find the lowest cherry.
            current_heights = copy_of_inputs.Update_Heights(
                current_heights, reducible_pairs
            )
            lowest_cherry = None
            lowest_height = None
            lowest_height_tuple = None
            lowest_heights_found = 1
            for pair in reducible_pairs:
                height_pair_tuple = current_heights[pair][0]
                height_pair = float(height_pair_tuple[0] + height_pair_tuple[1]) / 2
                new_found = False
                if (not lowest_height) or lowest_height > height_pair:
                    new_found = True
                    lowest_heights_found = 1
                elif lowest_height == height_pair:
                    lowest_heights_found += 1
                    if random.random() < 1 / float(lowest_heights_found):
                        new_found = True
                if new_found:
                    lowest_cherry = pair
                    lowest_height = height_pair
                    lowest_height_tuple = height_pair_tuple

            CPS += [lowest_cherry]
            heights_seq += [lowest_height_tuple]
            reduced_by_lowest_cherry = copy_of_inputs.Reduce_Pair_In_All(
                lowest_cherry, reducible_pairs=reducible_pairs
            )
            reducible_pairs = copy_of_inputs.Update_Reducible_Pairs(
                reducible_pairs, reduced_by_lowest_cherry
            )
            reduced_trees += [reduced_by_lowest_cherry]
            candidate_leaves = set(lowest_cherry)
        return CPS, reduced_trees, heights_seq

    def Update_Heights(self, current_heights, reducible_pairs):
        """
        Returns an updated dictionary of heights of the reducible pairs
        """
        for pair, trees in reducible_pairs.items():
            # updating is only necessary when the set of trees for that pair is changed or the reducible pair was not reducible before.
            if not pair in current_heights or not current_heights[pair][1] == len(
                trees
            ):
                height_pair = self.Height_Pair(pair, trees)
                current_heights[pair] = (height_pair, len(trees))
        return current_heights

    def Height_Pair(self, pair, trees):
        """
        Returns the average height of a pair in a set of trees
        The pair must be reducible in each tree in 'trees'
        """
        height_pair = [0, 0]
        for t in trees:
            height_in_t = cherry_height(self.trees[t], *pair)
            height_pair[0] += height_in_t[0]
            height_pair[1] += height_in_t[1]
        return [height_pair[0] / float(len(trees)), height_pair[1] / float(len(trees))]

    def Find_All_Pairs(self):
        """
        Finds the set of reducible pairs in all trees
        Returns a dictionary with reducible pairs as keys, and the trees they reduce as values.
        """
        reducible_pairs = dict()
        for i, t in self.trees.items():
            red_pairs_t = find_all_reducible_pairs(t)
            for pair in red_pairs_t:
                if pair not in reducible_pairs:
                    reducible_pairs[pair] = set()
                reducible_pairs[pair].add(i)
        return reducible_pairs

    def Update_Reducible_Pairs(self, reducible_pairs, new_red_trees):
        """
        Returns the updated dictionary of reducible pairs in all trees after a reduction (with the trees they reduce as values)
        we only need to update for the trees that got reduced: 'new_red_treed'
        """
        # Remove trees to update from all pairs
        for pair, trees in list(reducible_pairs.items()):
            trees.difference_update(new_red_trees)
            if len(trees) == 0:
                del reducible_pairs[pair]
        # Add the trees to the right pairs again
        for index in new_red_trees:
            if index in self.trees:
                t = self.trees[index]
                red_pairs_t = find_all_reducible_pairs(t)
                for pair in red_pairs_t:
                    if pair in reducible_pairs:
                        reducible_pairs[pair].add(index)
                    else:
                        reducible_pairs[pair] = set([index])
        return reducible_pairs

    def Reduce_Pair_In_All(self, pair, reducible_pairs=dict()):
        """
        Reduces the given pair in all networks in the problem.
        Returns the set of networks that were reduced.

        :note: This method changes the set of networks in the problem, so it should only be called in a copy of the problem.

        :param pair: a pair of leaves.
        :param reducible_pairs: a dictionary of reducible pairs and the networks in which they are reducible, as returned by Find_All_Pairs.
        """

        reduced_trees_for_pair = []
        if pair in reducible_pairs:
            trees_to_reduce = reducible_pairs[pair]
        else:
            # pair not found, trying all trees
            trees_to_reduce = deepcopy(self.trees)
        for i in trees_to_reduce:
            if i in self.trees:
                t = self.trees[i]
                t, cherry_type = reduce_pair(t, *pair, inplace=True)
                if cherry_type == CHERRYTYPE.RETICULATEDCHERRY:
                    reduced_trees_for_pair += [i]
                elif cherry_type == CHERRYTYPE.CHERRY:
                    reduced_trees_for_pair += [i]
                    t.leaves.remove(pair[0])
                if len(t.edges()) <= 1:
                    del self.trees[i]
        return set(reduced_trees_for_pair)

    def Reduce_Trivial_Pairs(self, candidate_leaves):
        """
        Reduces the trivial pairs in the current set of networks.
        Runs efficiently by giving a set of leaves 'candidate_leaves' that may be involved in trivial pairs.
        This set must be given; after a reduction of the pair (a,b) only using the leaves a and b works.
        Returns the reduced pairs and the sets of networks that were reduced.

        :note: This method changes the set of networks in the problem, so it should only be called in a copy of the problem.
        """
        seq = []
        reduced_tree_sets = []
        while candidate_leaves:
            l = candidate_leaves.pop()
            new_pairs = list(self.Trivial_Pair_With(l))
            if new_pairs:
                seq += new_pairs
                for p in new_pairs:
                    red_trees_p = self.Reduce_Pair_In_All(p)
                    reduced_tree_sets += [red_trees_p]
                    candidate_leaves = candidate_leaves | set(p)
        return seq, reduced_tree_sets

    def Reduce_Trivial_Pairs_Store_Pairs(self, candidate_leaves, reducible_pairs):
        """
        Reduces the trivial pairs in the current set of networks.
        Runs efficiently by giving a set of leaves 'candidate_leaves' that may be involved in trivial pairs.
        This set must be given; after a reduction of the pair (a,b) only using the leaves a and b works.
        Returns the reduced pairs and the sets of networks that were reduced, also updates the reducible pairs.

        This method is similar to Reduce_Trivial_Pairs, but it uses and updates the dictionary of reducible pairs to reduce the trivial pairs.

        :note: This method changes the set of networks in the problem, so it should only be called in a copy of the problem.
        """
        seq = []
        reduced_tree_sets = []
        while candidate_leaves:
            l = candidate_leaves.pop()
            new_pairs = list(self.Trivial_Pair_With(l))
            if new_pairs:
                seq += new_pairs
                for p in new_pairs:
                    red_trees_p = self.Reduce_Pair_In_All(
                        p, reducible_pairs=reducible_pairs
                    )
                    reducible_pairs = self.Update_Reducible_Pairs(
                        reducible_pairs, red_trees_p
                    )
                    reduced_tree_sets += [red_trees_p]
                    candidate_leaves = candidate_leaves | set(p)
        return seq, reduced_tree_sets, reducible_pairs

    def Reduce_Trivial_Pairs_Lengths(self, candidate_leaves, reducible_pairs):
        """
        Reduces the trivial pairs in the current set of networks.
        Runs efficiently by giving a set of leaves 'candidate_leaves' that may be involved in trivial pairs.
        This set must be given; after a reduction of the pair (a,b) only using the leaves a and b works.
        Returns the reduced pairs and the sets of networks that were reduced, also updates the reducible pairs and their heights.

        This method is similar to Reduce_Trivial_Pairs, but it uses and updates the dictionary of reducible pairs to reduce the trivial pairs.
        It also keeps track of the heights of the reduced pairs.

        :note: This method changes the set of networks in the problem, so it should only be called in a copy of the problem.
        """
        seq = []
        reduced_tree_sets = []
        heights_seq = []
        while candidate_leaves:
            l = candidate_leaves.pop()
            new_pairs = list(self.Trivial_Pair_With(l))
            if new_pairs:
                seq += new_pairs
                for p in new_pairs:
                    height_p = self.Height_Pair(p, reducible_pairs[p])
                    red_trees_p = self.Reduce_Pair_In_All(
                        p, reducible_pairs=reducible_pairs
                    )
                    heights_seq += [height_p]
                    reducible_pairs = self.Update_Reducible_Pairs(
                        reducible_pairs, red_trees_p
                    )
                    reduced_tree_sets += [red_trees_p]
                    candidate_leaves = candidate_leaves | set(p)
        return seq, reduced_tree_sets, reducible_pairs, heights_seq

    def Trivial_Pair_With(self, l):
        """
        Returns all trivial pairs involving the leaf l as first element of the pair.
        """
        pairs = set()
        # Go through all trees t with index i.
        for i, t in self.trees.items():
            # If the leaf occurs in t
            if l in t.leaves:
                # Compute reducible pairs of t with the leaf as first coordinate
                pairs_in_t = set(find_reducible_pairs_with_first(t, l))
                # If we did not have a set of candidate pairs yet, use pairs_in_t
                if not pairs:
                    pairs = pairs_in_t
                # Else, the candidate pairs must also be in t, so take intersection
                else:
                    pairs = pairs & set(pairs_in_t)
                # If we do not have any candidate pairs after checking a tree with l as leaf, we stop.
                if not pairs:
                    break
        return pairs

    def Improve_Sequence(self, CPS, reduced_trees, progress=False):
        """
        Improves a sequence 'CPS' for the input trees by removing elements and checking whether the new sequence still reduces all trees.
        Returns this improved sequence and the corresponding sets of reduced trees for each pair.
        """
        seq = deepcopy(CPS)
        i = 0
        while i < len(seq):
            redundant = True
            relevant_tree_indices = reduced_trees[i]
            new_relevant_pairs_for_trees = dict()
            for j in relevant_tree_indices:
                # Check if the shorter sequence reduces the trees, and if so, record which pairs reduced a cherry in which tree
                new_relevant_pairs_for_trees[j] = get_indices_of_reducing_pairs(
                    seq[:i] + seq[i + 1 :],
                    self.trees[j],
                )
                if not new_relevant_pairs_for_trees[j]:
                    redundant = False
                    break
            if redundant:
                # Remove the ith element from seq and reduced trees
                seq.pop(i)
                reduced_trees.pop(i)
                # Update reduced_trees for the relevant trees
                # First remove all j in relevant_tree_indices
                for tree_index_set in reduced_trees:
                    for j in relevant_tree_indices:
                        tree_index_set.discard(j)
                # Now add them back at the right places, according to "new_relevant_pairs_for_trees"
                for j in relevant_tree_indices:
                    for index in new_relevant_pairs_for_trees[j]:
                        reduced_trees[index].add(j)
                if progress:
                    print("New length is " + str(len(seq)))
                    print("Continue at position " + str(i))
            else:
                i += 1
        return seq, reduced_trees
