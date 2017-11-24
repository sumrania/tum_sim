
import random

import gene_prof as gnpr
import chrm_prof as chpr

mut_size_mean = 5 # mean mutation size

len_chrm_1 = len("AAABBBCCCDDDEEEFFF")
len_chrm_2 = len("RRRRRSSSSSTTTTTUUUUUVVVVV")
len_chrm_3 = len("XXXXYYYYZZZZ")

# build a dictionary of chromosome profiles before creating a gene profile.
#   key is tuple (s, i) where s (str) is chromosome name and i (int) is number of homologous chrm
#                example: ('1', 0) is maternal chromosome 1,     example: ('2', 1) is paternal chromosome 2
#   val is ChrmProf object
chrm_dict = {}
chrm_dict[('1', 0)] = chpr.ChrmProf(len_chrm_1)
chrm_dict[('1', 1)] = chpr.ChrmProf(len_chrm_1)

chrm_dict[('2', 0)] = chpr.ChrmProf(len_chrm_2)
chrm_dict[('2', 1)] = chpr.ChrmProf(len_chrm_2)

chrm_dict[('3', 0)] = chpr.ChrmProf(len_chrm_3)
chrm_dict[('3', 1)] = chpr.ChrmProf(len_chrm_3)

gp = gnpr.GeneProf(chrm_dict, mut_size_mean)
gp.mutate()
gp.pprint()
print ''

gp.mutate()
print '\n\n'
print 'gp info:'
gp.pprint()

print '\n\n'
gp_copied = gp.deepcopy()

gp_copied.mutate()
print 'gp copied info:'
gp_copied.pprint()

print '\n\n'
print 'gp info:'
gp.pprint()


