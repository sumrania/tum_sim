
import chrm_prof as chpr
import random
import numpy as np
import sys
import os


class GeneProf:

	# input: chrm_dict (dict) key is tuple (s, i) where s (str) is chromosome name and i (int) is number of homologous chrm
	#                           example: ('1', 0) is maternal chromosome 1,     example: ('2', 1) is paternal chromosome 2
	#                         val is ChrmProf object
	#        mut_size_mean (float) mean size of a random mutation
	#        mut_size_var (float) variance in size of random mutation
	def __init__(self, chrm_dict, mut_size_mean):
		self.mut_types = ['amp', 'rem', 'inv']
		self.chrm_dict = chrm_dict
		self.mut_size_mean = mut_size_mean

	# randomly mutates one chromosome with one random mutation
	def mutate(self):
		mut_type, chrm_tup, bgn, end = self._get_legal_random_mutation()
		if mut_type == 'amp':
			self.chrm_dict[chrm_tup].amp(bgn, end)
		elif mut_type == 'rem':
			self.chrm_dict[chrm_tup].rem(bgn, end)
		else:
			self.chrm_dict[chrm_tup].inv(bgn, end)

	# output: copy_num_dict (dict) key is chromosome (str). val is tuple of 3 lists
	#                key: chromosome (str)
	#                val: tuple: (bgns, ends, cps) of type (list of int, list of int, list of int)
	#                             bgns (list of int) beginning positions for each segment
	#                             ends (list of int) ending    positions for each segment
	#                             cps  (list of float) mixed copy number for each segment
	def get_copy_nums_dict(self):
		result = {}
		usages = [1, 1]
		for (chrm_str, pm) in self.chrm_dict:
			if chrm_str not in result:
				(bgns_m, ends_m, cps_m) = self.chrm_dict[(chrm_str, 0)].get_copy_nums()
				(bgns_p, ends_p, cps_p) = self.chrm_dict[(chrm_str, 1)].get_copy_nums()
				triplets = [[bgns_p, ends_p, cps_p], [bgns_m, ends_m, cps_m]]
				[res_bgns, res_ends, res_cps] = combine_copy_nums(triplets, usages)
				result[chrm_str] = (res_bgns, res_ends, res_cps)
		return result

	# deep copies GeneProf object and returns copy
	def deepcopy(self):
		chrm_dict_new = {}
		for (chrm_str, pm) in list(self.chrm_dict.keys()):
			chrm_dict_new[(chrm_str, pm)] = self.chrm_dict[(chrm_str, pm)].deepcopy()
		return GeneProf(chrm_dict_new, self.mut_size_mean)

	def pprint(self):
		l = sorted(self.chrm_dict.keys())
		for chrm_str, pm in l:
			_printnow('\n(' + str(chrm_str) + ', ' + str(pm) + ')\n')
			self.chrm_dict[chrm_str, pm].pprint()
		print 'copy_num_dict:', self.get_copy_nums_dict()

	#
	#   PRIVATE
	#

	# output: mut_type (str) one of 'amp', 'rem', 'inv'
	#         chrm_tup (tuple) (s, i) where s (str) is chromosome name and i (int) is number of homologous chrm
	#         bgn (int), end (int) are beginning and end positions on chromosome
	def _get_random_mutation(self):
		mut_type = random.choice(self.mut_types)
		chrm_tup = random.choice(list(self.chrm_dict.keys()))
		chrm_len = self.chrm_dict[chrm_tup].n
		
		mut_size = int(round(np.random.exponential(self.mut_size_mean)))
		if mut_size <= 0:         # mutation size must be positive number
			mut_size = 1
		elif mut_size > chrm_len: # mutation size cannot be larger than the chromosome
			mut_size = chrm_len

		bgn = random.randint(0, chrm_len - mut_size)
		end = bgn + mut_size - 1

		return mut_type, chrm_tup, bgn, end


	def _is_legal_mutation(self, chrm_tup, bgn, end):
		chrm_str = chrm_tup[0]
		chrm_tup_mat = (chrm_str, 0) # maternal chromosome tuple
		chrm_tup_pat = (chrm_str, 1) # paternal ...
		if bgn > end:
			return False
		if not self.chrm_dict[chrm_tup_mat]._is_splitable(bgn, end):
			return False
		if not self.chrm_dict[chrm_tup_pat]._is_splitable(bgn, end):
			return False
		return True


	# geneProf_list contains list of geneProfs
	def _get_legal_random_mutation(self):
		mut_type, chrm_tup, bgn, end = self._get_random_mutation()
		while not self._is_legal_mutation(chrm_tup, bgn, end):
			mut_type, chrm_tup, bgn, end = self._get_random_mutation()
		return mut_type, chrm_tup, bgn, end

	def print_chrm_seq(self):
		l = sorted(self.chrm_dict.keys())
		for (idx,pm) in l:
			print '(', idx, ',', pm, '): ', self.chrm_dict[(idx,pm)].chrm


	# input: cov (int), read_len (int)
	# output: sv_dict
	#         key: chromsome index (str)
	#         val: dictionary
	#              key: bp tuple (pos(int), isLeft (bool))
	#              val: dictionary
	#                   key: "mate", val: mate_bp tuple (pos(int), isLeft (bool))
	#                   key: "mated_reads", val: number of reads containing both curr and mated bp (int)
	#                   key: "total_reads", val: number of total reads (int)
	#                   key: "copy_num", val: copy number of current bp (int)
	def get_sv_read_nums_dict(self, cov, read_len):
		result = dict()
		for (idx, pm) in self.chrm_dict:
			if idx not in result:
				temp = dict()
				# paternal chrom
				sv_dict_p = self.chrm_dict[(idx, 0)].get_sv_read_nums(cov, read_len)

				# maternal chrom
				sv_dict_m = self.chrm_dict[(idx, 1)].get_sv_read_nums(cov, read_len)

				# combine paternal and maternal chrom sv dict
				repeated = set()
				for (pos, isLeft) in sv_dict_p:
					if (pos, isLeft) not in sv_dict_m:
						temp[(pos, isLeft)] = sv_dict_p[(pos, isLeft)]
					else:
						repeated.add((pos, isLeft))
						if sv_dict_p[(pos, isLeft)]["mate"] == sv_dict_m[(pos, isLeft)]["mate"]:
							temp[(pos, isLeft)]["mate"] = sv_dict_p[(pos, isLeft)]["mate"]
							temp[(pos, isLeft)]["mated_reads"] = sv_dict_p[(pos, isLeft)]["mated_reads"] + sv_dict_m[(pos, isLeft)]["mated_reads"]
							temp[(pos, isLeft)]["total_reads"] = sv_dict_p[(pos, isLeft)]["total_reads"] + sv_dict_m[(pos, isLeft)]["total_reads"]
							temp[(pos, isLeft)]["copy_num"] = sv_dict_p[(pos, isLeft)]["copy_num"] * 0.5 + sv_dict_m[(pos, isLeft)]["copy_num"] * 0.5
						else:
							print "\n"
							print 'Different mate bp in pair of chromosomes!!!'
							print (pos, isLeft), sv_dict_p[(pos, isLeft)]["mate"], sv_dict_m[(pos, isLeft)]["mate"]
							print "\n"
				for (pos, isLeft) in sv_dict_m:
					if (pos, isLeft) not in repeated:
						temp[(pos, isLeft)] = sv_dict_m[(pos, isLeft)]
				result[idx] = temp
		return result


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   C O M B I N E   C O P Y   N U M B E R   F U N C T I O N S   #
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

def _printnow(s):
	s = str(s)
	sys.stdout.write(s)
	sys.stdout.flush()

# type(triplets) = list, each triplet has 3 lists (bgns, ends, cps)
# type(usages) = list
#  input: triplets (3D list) [num_chromosomes by 3 by num_segments]
#                            triplets[k, 0, s] is beginning position for segment s on chromosome k
#                            triplets[k, 1, s] is ending    position for segment s on chromosome k
#                            triplets[k, 2, s] is copy number        for segment s on chromosome k
#         usages (list of float) [num_chromosomes] long. usages[k] is weight (percent) of chromosome k
# output: bgns (list of int) [num_unioned_segments] bgns[s] is beginning position of segment s
#         ends (list of int) [num_unioned_segments] ends[s] is ending    position of segment s
#         cps  (list of int) [num_unioned_segments] cps[s]  is fractional copy number of segment s
def combine_copy_nums(triplets, usages):
	posList, pos_dir_dict = _get_pos_info_dict(triplets)

	combined_segment_pos_list = _get_combined_segment_pos(posList, pos_dir_dict) # list of (bgnPos,endPos)

	bgn2idx, end2idx = _get_indices_dict(combined_segment_pos_list) # key: bgn or end position, val: idx of segment in combined chrom

	res_bgns, res_ends = list(), list()
	for (bgnPos,endPos) in combined_segment_pos_list:
		res_bgns.append(bgnPos)
		res_ends.append(endPos)

	res_cps = [0] * len(res_bgns)
	for i in range(len(triplets)):
		[bgns, ends, cps] = triplets[i]
		for j in range(len(bgns)):
			idx_list = _get_indices_for_segment(bgn2idx, end2idx, bgns[j], ends[j])
			for k in idx_list:
				res_cps[k] += float("{0:.2f}".format(cps[j] * usages[i]))

	return [res_bgns, res_ends, res_cps]
	

def _get_pos_info_dict(triplets):
	posSet = set() # combine bgn and end positions
	pos_dir_dict = dict() # key: pos, val: ['bgn'] or ['end'] or ['end','bgn']
	for i in range(len(triplets)):
		[bgns, ends, cps] = triplets[i]
		for j in range(len(bgns)):
			posSet.add(bgns[j])
			posSet.add(ends[j])
			if bgns[j] not in pos_dir_dict:
				pos_dir_dict[bgns[j]] = ['bgn']
			else:
				if 'bgn' not in pos_dir_dict[bgns[j]]:
					pos_dir_dict[bgns[j]] += ['bgn']
			if ends[j] not in pos_dir_dict:
				pos_dir_dict[ends[j]] = ['end']
			else:
				if 'end' not in pos_dir_dict[ends[j]]:
					pos_dir_dict[ends[j]] += ['end']
	posList = sorted(list(posSet))
	return posList, pos_dir_dict

# return (bgn, end) position for combined segment, eg. [(bgn1, end1), (bgn2, end2), ...]
def _get_combined_segment_pos(posList, pos_dir_dict):
	result = list()
	tempBgn = posList[0]
	i = 0

	while i < len(posList) - 1:
		# end at this position
		if 'end' in pos_dir_dict[posList[i]]:
			tempEnd = posList[i]
			result.append((tempBgn, tempEnd))

			if 'bgn' in pos_dir_dict[posList[i + 1]]:
				tempBgn = posList[i + 1]
			else:
				tempBgn = posList[i] + 1
			i += 1
		# not end at this position
		else:
			if 'bgn' in pos_dir_dict[posList[i + 1]] and 'end' not in pos_dir_dict[posList[i + 1]]:
				tempEnd = posList[i + 1] - 1
				result.append((tempBgn, tempEnd))
				tempBgn = posList[i + 1]
				i += 1
			elif 'end' in pos_dir_dict[posList[i + 1]] and 'bgn' not in pos_dir_dict[posList[i + 1]]:
				tempEnd = posList[i + 1]
				result.append((tempBgn, tempEnd))
				if tempEnd == posList[-1]:
					break
				else:
					tempBgn = posList[i + 2]
				i += 2
			else:
				tempEnd = posList[i + 1] - 1
				result += [(tempBgn, tempEnd), (posList[i + 1], posList[i + 1])]
				tempBgn = posList[i + 1] + 1
				i += 1

	# last element in posList
	if tempBgn == posList[-1]:
		result.append((tempBgn, tempBgn))
	return result


def _get_indices_dict(combined_segment_pos_list):
	bgn2idx = dict()
	end2idx = dict()
	idx = 0
	for (s,e) in combined_segment_pos_list:
		bgn2idx[s] = idx
		end2idx[e] = idx
		idx += 1
	return bgn2idx, end2idx


# given start and end position, output list of segment indices (continuous)
def _get_indices_for_segment(bgn2idx, end2idx, s, e):
	result = list()
	firstIdx = bgn2idx[s]
	lastIdx = end2idx[e]
	for i in range(firstIdx, lastIdx + 1):
		result.append(i)
	return result
