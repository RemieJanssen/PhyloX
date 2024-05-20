.. PhyloX documentation master file, created by
   sphinx-quickstart on Sun Jul  9 20:29:15 2023.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to PhyloX's documentation!
==================================

Phylox DiNetwork
----------------

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:
   :caption: Phylox DiNetwork

   phylox.dinetwork
   phylox.newick_parser
   phylox.isomorphism

Network properties
------------------

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:
   :caption: Network properties

   phylox.classes
   phylox.networkproperties

Network operations
------------------

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:
   :caption: Network operations

   phylox.generators
   phylox.cherrypicking
   phylox.rearrangement

Supporting modules
------------------

.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:
   :caption: Supporting modules

   phylox.constants
   phylox.exceptions


PhyloX is a python package for parsing, manipulating, and analysing phylogenetic networks.
By building upon the widely used `NetworkX <https://networkx.github.io/>`_ package, PhyloX provides a simple and intuitive interface for working with phylogenetic networks.
Unlike other packages like `Biopython <https://biopython.org/>`_ and `DendroPy <https://dendropy.org/>`_, PhyloX is designed to work with phylogenetic networks, rather than phylogenetic trees.
Moreover, PhyloX is designed to work with networks that have internally labelled nodes, multiple roots, and non-binary nodes (although not all methods are implemented for all these types of networks).

Some of the features of PhyloX include:
 - Generating or reading and writing phylogenetic networks.

   - Reading from extended Newick format (:func:`phylox.dinetwork.DiNetwork.from_newick`) with internal node labels, edge lengths, support values, and inheritance probabilities.
   - Building from a cherry-picking sequence (:func:`phylox.dinetwork.DiNetwork.from_cherry_picking_sequence`)
   - Generating random phylogenetic networks with several different models (:mod:`phylox.generators`).

 - Modifying and analysing phylogenetic networks using rearrangement moves (:mod:`phylox.rearrangement`).

   - Calculating the rearrangement distance between two networks (exact and heuristically).
   - Generating the set of all networks that can be reached from a given network using a single rearrangement move.
   - Moves are robust objects, that give meaningful error messages when they are applied to networks that they are not applicable to.

 - Cherry-picking methods for phylogenetic networks (:mod:`phylox.cherrypicking`).

   - Checking whether a network is an orchard network (:func:`phylox.classes.is_orchard`). 
   - Network containment checking for tree-child networks using cherry-picking sequences.
   - Combining networks using cherry-picking sequences.

 - Isomorpism functions for phylogenetic networks (:mod:`phylox.isomorphism`).

   - Checking whether two networks are isomorphic.
   - Counting the number of automorphisms of a network.

A use-case that neatly combines many of these features is the following. Suppose we want a test set of orchard networks with 10 leaves and 5 reticulations. These networks can be generated randomly with a Metropolis-Hasting sampling.

We start with an arbitrary orchard network with 10 leaves and 5 reticulations, which can be generated with the function :func:`phylox.generators.randomTC:generate_network_random_tree_child_sequence`.
Then repeatedly apply some (large number of) random rearrangement moves, only actually applying moves if the resulting network is orchard.
This can be done with the built-in MCMC sampler, :func:`phylox.generators.mcmc:sample_mcmc_networks`.
To ensure all sampled networks are orchard, we use the `restriction_map` argument to specify that only orchard networks should be sampled, and to ensure the right number of leaves and reticulations are sampled, we set the `move_type_probabilities` for vertical moves to 0.
We can use the `sample_size` argument to specify how many networks we want to sample.
Finally, we write the sampled networks to a file in newick format.

The final code looks like this:

.. code-block:: python

   from phylox.generators.randomTC import generate_network_random_tree_child_sequence
   from phylox.generators.mcmc import sample_mcmc_networks
   from phylox.classes import is_orchard
   from phylox.rearrangement.move import MoveType

   # Generate an arbitrary orchard network with 10 leaves and 5 reticulations
   start_network = generate_network_random_tree_child_sequence(10, 5, seed=4321)
   # Generate 100 orchard networks with 10 leaves and 5 reticulations
   sampled_networks = sample_mcmc_networks(
      start_network,
      {MoveType.TAIL: 0.5, MoveType.HEAD: 0.5},
      number_of_samples=100,
      burn_in=5,
      restriction_map=is_orchard,
      seed=1234,
      add_root_if_necessary=True,
      correct_symmetries=False,
   )
   # Write the sampled networks to a file
   with open("sampled_networks.nwk", "w") as f:
      for network in sampled_networks:
         f.write(network.newick() + "\n")


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
