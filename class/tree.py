
import sys   # hack for printing special characters
reload(sys)
sys.setdefaultencoding('utf-8')

import chrm_prof as cp
import gene_prof as gp
import numpy as np
import math

from graphviz import Digraph
from colour import Color # for coloring nodes in trees

import treelib

class Tree:
	# chrm_lens (list of int) length of each chromosome
	# mut_size_mean (float) mean size of random mutation
	def __init__(self, chrm_lens, mut_size_mean, alpha, beta):
		self.root = Node(gp.GeneProf(_get_chrm_dict(chrm_lens), mut_size_mean), 1.0, None, 0.0, 0)
		self.leaves = [self.root]
		self.time = 0.0
		self.alpha = alpha
		self.beta = beta

	def mutate(self):
		rates = [ leaf.pop for leaf in self.leaves ]
		times = [ np.random.exponential(1.0 / rate) for rate in rates ]
		leaf_i = np.argmin(times)
		leaf = self.leaves[leaf_i]
		self.time += times[leaf_i]

		# create children of leaf to be mutated
		gene_prof = leaf.gene_prof
		gene_prof_left = gene_prof.deepcopy()
		gene_prof_left.mutate()
		gene_prof_right = gene_prof.deepcopy()

		# distribute percent of cell types to children
		percent_left = np.random.beta(self.alpha, self.beta)
		leaf.left = Node(gene_prof_left, leaf.pop * percent_left, leaf, self.time, leaf.num_muts + 1)
		leaf.right = Node(gene_prof_right, leaf.pop * (1 - percent_left), leaf, leaf.time, leaf.num_muts)

		# record new children as leaves. remove old leaf from list of leaves
		self.leaves.pop(leaf_i)
		self.leaves.append(leaf.left)
		self.leaves.append(leaf.right)

	# returns usage (list of float) random generated percent of each cell type
	def get_rand_sample(self, num_trials):
		if num_trials < 1:
			num_trials = 1
		pops = [ leaf.pop for leaf in self.leaves ]
		leaf_counts = np.zeros(len(self.leaves), dtype = float)
		for _ in xrange(0, num_trials):
			leaf = np.random.choice(self.leaves, p = pops)
			i = self.leaves.index(leaf)
			leaf_counts[i] += 1
		return list(leaf_counts / round(np.sum(leaf_counts)))

	def print_img(self, dname):
		bgnColor = Color("#fff200")
		colors = list(bgnColor.range_to(Color("#ff3d3d"), 1 * len(self.leaves)))
		tot_time = self.time

		self._label()
		dot = Digraph(format = 'png')
		stack = [self.root]
		while stack:
			cur = stack.pop()

			# calculate color from time created
			idx = int(math.floor(float(len(colors)) * cur.time / tot_time))
			if idx >= len(colors):
				idx = len(colors) - 1
			color = colors[idx]

			dot.attr('node', style='filled', color = color.hex_l)
			dot.node(str(id(cur)), label = str(cur._calc_L1_dist()))
			if cur.left:
				stack.append(cur.left)
				dot.edge(str(id(cur)), str(id(cur.left)), label = '')
			if cur.right:
				stack.append(cur.right)
				dot.edge(str(id(cur)), str(id(cur.right)), label = '')
		dot.render(dname + 'expanded_tree')

		# collapsed tree
		dot = Digraph(format = 'png')
		stack = [self.root]
		leaf_labels_seen = []
		while stack:
			cur = stack.pop()
			# calculate color from time created
			idx = int(math.floor(float(len(colors)) * cur.time / tot_time))
			if idx >= len(colors):
				idx = len(colors) - 1
			color = colors[idx]
			dot.attr('node', style='filled', color = color.hex_l)
			if cur.label not in leaf_labels_seen:
				leaf_labels_seen.append(cur.label)
				dot.node(cur.label, label = str(cur._calc_L1_dist()))
			if cur.left:
				stack.append(cur.left)
				dot.edge(cur.label, cur.left.label)
			if cur.right:
				stack.append(cur.right)
				# dot.edge(cur.label, cur.right.label)
		dot.render(dname + 'tree')

	def pprint(self):
		self._label()
		my_treelib_tree = treelib.Tree()
		queue = [self.root]
		while queue:
			cur = queue.pop()
			if cur.left:
				queue.append(cur.left)
			if cur.right:
				queue.append(cur.right)
			if cur.parent:
				my_treelib_tree.create_node(str(cur), str(cur.pop), parent = str(cur.parent.pop))
			else:
				my_treelib_tree.create_node(str(cur), str(cur.pop))
		my_treelib_tree.show()

	# labels all nodes from leaves
	def _label(self):
		for i, leaf in enumerate(self._get_leaves_right_to_left()):
			leaf.label = str(i)
			cur = leaf
			while cur.parent and cur.parent.right == cur:
				cur = cur.parent
				cur.label = cur.right.label

	def _get_leaves_right_to_left(self):
		leaves = []
		stack = [self.root]
		while stack:
			cur = stack.pop()
			if cur.left:
				stack.append(cur.left)
			if cur.right:
				stack.append(cur.right)
			if not cur.left and not cur.right: # at a leaf
				leaves.append(cur)
		return leaves

class Node:
	def __init__(self, gene_prof, pop, par, time_created, num_muts):
		self.gene_prof = gene_prof
		self.pop = pop
		self.left = None
		self.right = None
		self.parent = par
		self.label = ''
		self.time = time_created
		self.num_muts = num_muts

	def __str__(self):
		return str(self.label)

	# calculates the L1 distance of the node's gene profile to the normal profile: [2, 2, ..., 2]
	def _calc_L1_dist(self):
		dist = 0
		for _, (_, _, cps) in self.gene_prof.get_copy_nums_dict().iteritems():
			dist += sum([ abs(cp - 2) for cp in cps ])
		return dist

def _get_chrm_dict(chrm_lens):
	chrm_dict = {}
	for i, chrm_len in enumerate(chrm_lens):
		chrm_dict[(str(i), 0)] = cp.ChrmProf(chrm_len)
		chrm_dict[(str(i), 1)] = cp.ChrmProf(chrm_len)
	return chrm_dict
