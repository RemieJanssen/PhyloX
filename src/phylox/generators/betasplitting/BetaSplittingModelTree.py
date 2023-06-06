############################################
# Beta-splitting model (Aldous1996)
#
# To run this script, type in your terminal:
# python3 BetaSplittingModel.py
#
# Parameters can be adjusted inside the code (line 74).
#
############################################

import os
import sys
import math
import random
import numpy as np
from scipy.special import gamma, loggamma
import dendropy
from pathlib import Path
import multiprocessing
from multiprocessing import Manager

############################################
# Simulation functions
############################################

# a_n is a normalizing constant defined in 
# Equation (2) of Aldous1996 (so the sum of 
# the values is equal to 1. It is not 
# computed here to save time.
def a_n(beta):
	return 1

# Compute the "probability" to split n in (i|n-1), where i=1,..,n-1
def computeProb(n,beta):
	q_n = []
	for i in range(1,n):
		q_i_n = np.exp((loggamma(beta+i+1)+loggamma(beta+n-i+1))-((a_n(beta)+loggamma(i+1)+loggamma(n-i+1))))
		q_n.append(q_i_n)
	return q_n

# n: number of tips
def simulateBetaSplitting(n, beta):
	# Initialize tree.
	taxon_namespace = dendropy.TaxonNamespace(["N" + str(i) for i in range(1, 2*n)])
	tree            = dendropy.Tree(taxon_namespace=taxon_namespace) # dendropy.Tree()
	setattr(tree.seed_node, "n", n)
	queue = [tree.seed_node]
	# Insert one node at each iteration.
	while queue:
		node       = queue.pop(0)
		node.taxon = taxon_namespace.get_taxon("N" + str(len(tree.nodes())))		
		n          = node.n
		# Internal node. Splits again.
		if (n > 1): 
			# Compute the "probability" to split n in (i|n-1), where i=1,..,n-1
			q_n   = computeProb(n,beta)
			split = random.choices(population=list(range(1,n)), weights=q_n, k=1)[0]
			# Create children.
			ch1 = node.new_child(edge_length=1)
			ch2 = node.new_child(edge_length=1)
			setattr(ch1, "n", split)
			setattr(ch2, "n", n-split)
			# Add children to the queue.
			queue.append(ch1)
			queue.append(ch2)
	# Return tree.
	return tree

############################################
# MAIN
############################################

############################################
# General parameters.

# All desired tree sizes.
tree_size_all   = [500,750,1000] # Trees with 500, 750, and 1000 tips respectively.
# Number of simulated trees per tree size.
nb_trees        = 10
# Range of beta parameter.
min_beta        = -2.0
max_beta        = 10.0
# Flag to indicate if files should overwrite possible existent files.
saveifexists    = False
# Number of cores in your computer (trees are simulated in parallel to save time).
nb_threads      = 4

############################################
# Simulate trees

# Each tree size
for tree_size in tree_size_all:
	sim_prefix =  "betaSplitting-tips" + str(tree_size)
	if (not saveifexists) and any(sim_file.startswith(sim_prefix) for sim_file in os.listdir(".")):
		print("WARNING! Skipping simulation (files already exists): tree_size=" + str(tree_size))
	else:

		# Method to simulate a single tree.
		def runsubproc(pars):
			# Simulate tree.
			tree_size, beta = pars
			tree            = simulateBetaSplitting(tree_size, beta)
			infostr = str(tree_size) + "\t" + str(beta) + "\n"
			treeStr = tree.as_string(schema="newick")
			return (infostr, treeStr)

		# Simulate trees.
		current_points = 0
		while (current_points < nb_trees):

			# Uniform sample parameters (beta values).
			print("Current number of points: " + str(current_points) + " / Desired: " + str(nb_trees))
			points_i        = 0
			beta_values_all = []
			while points_i < (nb_trees-current_points):
				beta = np.random.uniform(low=min_beta, high=max_beta)
				beta_values_all.append(beta)
				points_i += 1

			# Simulate
			all_sim   = [(tree_size, beta) for beta in beta_values_all]

			print("  Run all simulations (" + str(len(all_sim)) +  ")...")
			p = multiprocessing.Pool(nb_threads)
			all_trees = p.map(runsubproc, all_sim)

			filepath_trees = sim_prefix + ".nwk"
			filepath_pars  = sim_prefix + ".tsv"

			if os.path.exists(filepath_trees) and os.path.exists(filepath_pars):
				append_write = 'a' # append if already exists
			else:
				append_write = 'w' # append if already exists
			
			with open(filepath_trees, append_write) as f:
				for (infostr, treeStr) in all_trees:
					if treeStr:
						f.write(treeStr)
						current_points += 1
				print("\t Saved: '" + str(filepath_trees) + "'")

			with open(filepath_pars, append_write) as f:
				for (infostr, treeStr) in all_trees:
					if treeStr:
						f.write(infostr)
				print("\t Saved: '" + str(filepath_pars) + "'")

