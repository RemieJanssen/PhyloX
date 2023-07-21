"""
A module for rearranging phylogenetic networks using rearrangement moves.

This module contains the class RearrangementProblem which is used to represent a rearrangement problem.
That is, a problem of transforming one phylogenetic network into another using rearrangement moves.
A solution to a rearrangement problem is a sequence of moves, which is a list of phylox.rearrangement.move.Move.

There are different types of horizontal and vertical moves, enumerated in phylox.rearrangement.move.MoveType.
For an overview of types of moves, and some results about the existence of solutions to rearrangement problems,
see for example "Rearranging Phylogentic Networks" by Remie Janssen.
The code for this module was originally written by Remie Janssen for the above mentioned thesis.
"""
