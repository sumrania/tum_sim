#     file: gen_data.py
#   author: Eaton Park, Murtazaflaza, Shefaroni
#  created: 12/4/2017
# modified: 12/4/2017
#  purpose: generates entire data set across many parameters by iteratively calling sim.py


# # # # # # # # # # #
#   I M P O R T S   #
# # # # # # # # # # #

import os
import sys
import argparse

# custom modules
import sim


# # # # # # # # # # # # #
#   C O N S T A N T S   #
# # # # # # # # # # # # #

DIRNAME = '/F/school/2018_carnegie_mellon/Semester_3_Fall_2017/Simulation\ of\ Bio/project/data/sim_data/'

DATA_SETS = {
	"toy": {
		"chrm_lens": [10000, 15000, 20000, 25000, 30000],
		"mean_mut_lens": [100, 1000]
	},
	"real": {
		"chrm_lens": [248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717, 133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285, 58617616, 64444167, 46709983, 50818468],
		"mean_mut_lens": [50000, 100000]
	}
} # real is 22 chromosomes excluding X and Y chromosomes
NUM_MUTS = [10, 50, 300]
NUM_SAMPS = [1, 3, 5, 10]
SAMP_VARS = {
	"low": 0.002,
	"med": 0.001,
	"high": 0.05
}
ALPHA_BETAS = [(5, 1), (1, 1)]



# # # # # # # # # # # # #
#   F U N C T I O N S   #
# # # # # # # # # # # # #

def main(argv):
	args = get_args(argv)
	dname = args['output_directory']
	print 'DO NOT RUN THIS UNLESS YOU WANT TO GENERATE AN ENTIRE DATA SET'
	print '\nREMOVE exit() CALL IN ACTUAL CODE TO CONTINUE'
	exit()

	for set_name, val in DATA_SETS.iteritems():
		chrm_lens = val['chrm_lens']
		mean_mut_lens = val['mean_mut_lens']
		for l in mean_mut_lens:
			for m in NUM_MUTS:
				for s in NUM_SAMPS:
					for vname, v in SAMP_VARS.iteritems():
						for a, b in ALPHA_BETAS:
							subdname = '_'.join([set_name, 'l', str(l), 'm', str(m), 's', str(s), 'v', vname, 'a', str(a), 'b', str(b)])
							tmp_dname = directorize(dname + subdname)
							make_dir(tmp_dname)
							printnow(tmp_dname)
							sim.sim(tmp_dname, chrm_lens, l, m, s, v, a, b)

def get_args(argv):
	parser = argparse.ArgumentParser(prog = 'gen_data.py', description = "generates entire data set across many parameters by iteratively calling sim.py")
	parser.add_argument('-o', '--output_directory', required = True, help = 'empty directory where all subdirectories will be created. each subdirectory will contain results of of a simulation run with different parameters', type = lambda x: valid_dir(parser, x))
	return vars(parser.parse_args(argv))

def valid_dir(parser, arg):
	if not os.path.exists(arg):
		parser.error('The directory \"' + str(arg) + '\" could not be found.')
	return directorize(arg)

def directorize(dname):
	if dname.endswith('/'):
		return dname
	return dname + '/'

def make_dir(dname):
	print dname
	if not os.path.exists(dname):
		os.makedirs(dname)

def printnow(s, newline = True):
	s = str(s)
	if newline:
		s += '\n'
	sys.stdout.write(s)
	sys.stdout.flush()


# # # # # # # # # # # # # # # # # # # # # # # # #
#   C A L L   T O   M A I N   F U N C T I O N   #
# # # # # # # # # # # # # # # # # # # # # # # # #

if __name__ == "__main__":
	main(sys.argv[1:])

