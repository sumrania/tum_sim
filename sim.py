#     file: sim.py
#   author: Eaton Park, Murtazaflaza, Shefaroni
#  created: 12/1/2017
# modified: 12/1/2017
#  purpose: generate tumor phylogeny using continuous time markov model for segmental mutations


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import sys      # for command line arguments
import os       # for manipulating files and folders
import argparse # for command line arguments
import numpy as np

# custom imports
sys.path.insert(0, 'class/')
sys.path.insert(0, 'help/')
import tree
import combine_segments as cs


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

# here


# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	sim(args['output_directory'], args['chrm_lens'], args['mean_mutation_length'], args['num_mutations'], args['num_samples'], args['sample_variance'], args['alpha'], args['beta'])

def sim(out_dir, chrm_lens, mean_mut_len, num_muts, num_samps, samp_var, alpha, beta):
	t = tree.Tree(chrm_lens, mean_mut_len, alpha, beta)
	for _ in xrange(0, num_muts):
		t.mutate()
	U = get_U(t, num_samps, samp_var)
	C = get_C(t)
	F = np.matmul(U, C)
	np.savetxt(out_dir + 'U.tsv', U, delimiter = '\t', fmt = '%1.6f')
	np.savetxt(out_dir + 'C.tsv', C, delimiter = '\t', fmt = '%1.6f')
	np.savetxt(out_dir + 'F.tsv', F, delimiter = '\t', fmt = '%1.6f')
	t.print_img(out_dir)

def get_U(t, num_samps, samp_var):
	num_leaves = len(t.leaves)
	num_muts = num_leaves - 1
	U = np.empty((num_samps, num_leaves), dtype = float)
	U[0, :] = [ leaf.pop for leaf in t._get_leaves_right_to_left() ] # first sample is perfect sample
	num_trials = int(round(float(num_muts) / ( samp_var * float(num_muts ** 2 + 2*num_muts + 1) )))
	for i in xrange(1, num_samps):
		U[i, :] = t.get_rand_sample(num_trials)
	return U

def get_C(t):
	num_leaves = len(t.leaves)

	X_bgns_dict = {} # key: chrm (str), val: X_bgns (2D list of int) X_bgns[c, s] is beginning position of segment s on chromosome c
	X_cps_dict = {}  # key: chrm (str), val: X_cps  (2D list of int) X_cps[c, s] is the copy number of segment s on chromosome c
	for leaf in t._get_leaves_right_to_left():
		cp_dict = leaf.gene_prof.get_copy_nums_dict()
		for chrm, (bgns, ends, cps) in cp_dict.iteritems():
			if chrm not in X_bgns_dict:
				X_bgns_dict[chrm] = []
				X_cps_dict[chrm] = []
			X_bgns_dict[chrm].append(bgns)
			X_cps_dict[chrm].append(cps)
	
	X_bgns_dict_out, X_cps_dict_out = {}, {}
	for chrm, X_bgns in X_bgns_dict.iteritems():
		X_cps = X_cps_dict[chrm]
		X_cps_dict_out[chrm], bgns = cs.comb_segs(X_bgns, X_cps)
		X_bgns_dict_out[chrm] = []
		for i in xrange(0, num_leaves):
			X_bgns_dict_out[chrm].append(bgns)

	num_segs = 0
	for _, X_bgns in X_bgns_dict_out.iteritems():
		num_segs += len(X_bgns[0])

	C = np.empty((num_leaves, num_segs), dtype = int)

	seg_indx = 0
	for chrm, X_cps in X_cps_dict_out.iteritems():
		chrm_num_segs = len(X_cps[0])
		for i in xrange(0, num_leaves):
			for j in xrange(seg_indx, seg_indx + chrm_num_segs):
				C[i, j] = X_cps[i][j - seg_indx]
		seg_indx += len(X_cps[0])
	
	return C


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M M A N D   L I N E   A R G U M E N T   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'sim.py', description = "generate tumor phylogeny using continuous time markov model for segmental mutations")
	parser.add_argument('chrm_lens', nargs = '*', help = 'space separated list of the lengths of each chromosome', type = lambda x: is_list_of_int(parser, x))
	parser.add_argument('-m', '--num_mutations', help = 'number of mutations that will occur in the tree', type = int, required = True)
	parser.add_argument('-l', '--mean_mutation_length', help = 'average size of a mutation', type = int, default = 500)
	parser.add_argument('-s', '--num_samples', help = 'number of tumor samples to take from tumor', type = int, default = 5)
	parser.add_argument('-v', '--sample_variance', help = 'variance in how cell types are chosen for sampling', type = float, default = 0.01)
	parser.add_argument('-a', '--alpha', type = lambda x: float_above(parser, x, 0.0), default = 5, help = 'alpha parmeter to beta distribution that determines the percent of cells that will belong to the mutated child for each mutation. smaller alpha increases population of normal child. alpha = beta = 1 is Uniform(0, 1)')
	parser.add_argument('-b', '--beta', type = lambda x: float_above(parser, x, 0.0), default = 2, help = 'beta parmeter to beta distribution that determines the percent of cells that will belong to the mutated child for each mutation. smaller beta increases population of mutant child. alpha = beta = 1 is Uniform(0, 1)')
	parser.add_argument('-o', '--output_directory', required = True, help = 'directory where all results of simulation will be placed', type = lambda x: valid_dir(parser, x))
	return vars(parser.parse_args(argv))

def is_list_of_int(parser, arg):
	try:
		arg = int(arg)
	except:
		parser.error(str(arg) + ' must be an integer')
	return arg

def float_above(parser, arg, lb):
	try:
		arg = float(arg)
	except:
		parser.error(str(arg) + ' must be a float')
	if arg <= lb:
		parser.error(str(arg) + ' must be strictly above ' + str(lb))
	return arg

def float_between(parser, arg, lb, ub):
	try:
		arg = float(arg)
	except:
		parser.error(str(arg) + ' must be a float')
	if arg < lb or arg > ub:
		parser.error(str(arg) + ' must be between ' + str(lb) + ' and ' + str(ub) + ' inclusive')
	return arg

def valid_dir(parser, arg):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	return directorize(arg)

def directorize(dname):
	if dname.endswith('/'):
		return dname
	return dname + '/'


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])
