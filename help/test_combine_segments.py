
import combine_segments as cs

X_bgns = []
X_bgns.append([1, 51, 76, 106])
X_bgns.append([1, 76, 151])
X_bgns.append([1, 66, 121])

X_cps = []
X_cps.append([2, 1, 3, 4])
X_cps.append([2, 3, 2])
X_cps.append([0, 2, 7])

Y_cps, bgns = cs.comb_segs(X_bgns, X_cps)
print bgns
print Y_cps