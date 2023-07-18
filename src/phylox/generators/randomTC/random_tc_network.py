import random

from phylox import DiNetwork

### For the tree child network containment paper by Murakami and Janssen
### This file contains functions to generate random tree-child (sub)networks


def random_tree_child_sequence(leaves, reticulations):
    """
    Returns a random tree-child sequence with a given number of leaves and reticulations
    :param leaves: number of leaves
    :param reticulations: number of reticulations
    :return: a random tree-child sequence

    :example:
    >>> from phylox.generators.randomTC.random_tc_network import random_tree_child_sequence
    >>> seq = random_tree_child_sequence(5, 2)
    >>> len(seq)
    6
    """
    current_leaves = set([1, 2])
    seq = [(2, 1)]
    not_forbidden = set([2])
    leaves_left = leaves - 2
    retics_left = reticulations

    # Continue until we added enough leaves and reticulations
    while leaves_left > 0 or retics_left > 0:
        # Decide if we add a leaf, or a reticulation
        type_added = "L"
        if len(not_forbidden) > 0 and leaves_left > 0 and retics_left > 0:
            if (
                random.randint(0, leaves_left + retics_left - 1) < retics_left
            ):  # probability of retic depends on number of reticulations left to add
                # if random.randint(0 , 1)<1:   #probability of reticulations and leaves are the same if both are an option
                type_added = "R"
        elif len(not_forbidden) > 0 and retics_left > 0:
            type_added = "R"
        elif leaves_left > 0:
            type_added = "L"
        else:
            return False

        # Actually add the pair
        if type_added == "R":
            first_element = random.choice(list(not_forbidden))
            retics_left -= 1
        if type_added == "L":
            first_element = len(current_leaves) + 1
            leaves_left -= 1
            current_leaves.add(first_element)
            not_forbidden.add(first_element)

        second_element = random.choice(list(current_leaves - set([first_element])))
        not_forbidden.discard(second_element)
        seq.append((first_element, second_element))

    # reverse the sequence, as it was built in reverse order
    seq = [pair for pair in reversed(seq)]
    return seq


def random_tree_child_subsequence(seq, r):
    """
    Returns a random tree-child subsequence with a given number of reticulations
    :param seq: a tree-child sequence
    :param r: number of reticulations in the subsequence
    :return: a random tree-child subsequence

    :example:
    >>> from phylox.generators.randomTC.random_tc_network import random_tree_child_subsequence
    >>> seq = [(4,1), (4, 1), (3, 1), (2, 1), (2, 1), (2,1)]
    >>> newSeq = random_tree_child_subsequence(seq, 1)
    >>> len(newSeq)
    4
    """
    # First `uniformly at random' choose one pair per leaf, with that leaf as first element
    leaves = dict()
    indices = set()
    for i, pair in enumerate(seq):
        x = pair[0]
        if x not in leaves:
            indices.add(i)
            leaves[x] = (1, i)
        else:
            if random.randint(0, leaves[x][0]) < 1:
                indices.remove(leaves[x][1])
                indices.add(i)
                leaves[x] = (leaves[x][0] + 1, i)
            else:
                leaves[x] = (leaves[x][0] + 1, leaves[x][1])
    # Add r reticulations with a max of the whole sequence
    unused = set(range(len(seq))) - indices
    for j in range(r):
        new = random.choice(list(unused))
        unused = unused - set([new])
        indices.add(new)
    newSeq = []
    for i, pair in enumerate(seq):
        if i in indices:
            newSeq.append(pair)
    return newSeq


def generate_network_random_tree_child_sequence(
    leaves, reticulations, label_leaves=False
):
    """
    Returns a random tree-child network with a given number of leaves and reticulations
    :param leaves: number of leaves
    :param reticulations: number of reticulations
    :return: a random tree-child network

    :example:
    >>> from phylox.generators.randomTC.random_tc_network import generate_network_random_tree_child_sequence
    >>> network = generate_network_random_tree_child_sequence(5, 2)
    >>> len(network.leaves)
    5
    >>> network.reticulation_number
    2
    """
    if leaves < 2 or reticulations < 0:
        raise ValueError(
            "Invalid number of leaves or reticulations, must be at least 2 and 0 respectively"
        )
    seq = random_tree_child_sequence(leaves, reticulations)
    if not seq:
        return False
    return DiNetwork.from_cherry_picking_sequence(seq, label_leaves=label_leaves)
