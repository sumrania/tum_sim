
# ALL CHROMOSOMES ASSUMED TO HAVE EQUAL LENGTH
#  input: X_bgns (2D list of int) X_bgns[c, s] is beginning position of segment s on chromosome c
#         X_cps  (2D list of int) X_cps[c, s] is the copy number of segment s on chromosome c
# output: out_X_cps (2D list of int) same as X_cps but now split into all segments from each chromosome
#         out_bgns (list of int) list of all common beginning positions
def comb_segs(X_bgns, X_cps):
	num_chrm = len(X_bgns)
	Xi = [ 0 for _ in X_bgns ] # keep track of where we are in each list
	out_bgns = []
	out_X_cps = [ [] for _ in X_bgns ]
	
	cur_bgns = ['Russell ROX #YangYang']
	while cur_bgns:
		cur_bgns = []
		for chrm, chrm_bgns in enumerate(X_bgns):
			cur_i = Xi[chrm]
			if cur_i < len(chrm_bgns):          # make sure we do not index X_bgns out of bounds on a chromosome
				cur_bgns.append(chrm_bgns[cur_i])

		if not cur_bgns:
			return out_X_cps, out_bgns

		min_bgn = min(cur_bgns)
		out_bgns.append(min_bgn)

		for chrm, chrm_bgns in enumerate(X_bgns):
			cur_i = Xi[chrm]
			if cur_i < len(chrm_bgns):      # make sure we do not index X_bgns out of bounds on a chromosome
				cur_bgn = chrm_bgns[cur_i]
				if min_bgn == cur_bgn:
					Xi[chrm] += 1
		for chrm, cps in enumerate(X_cps):
			out_X_cps[chrm].append(cps[Xi[chrm] - 1])
	return out_X_cps, out_bgns
