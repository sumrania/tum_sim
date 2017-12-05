import tree

t = tree.Tree([10000, 20000, 30000], 500, 1.0, 1.0)
# t.pprint()

for i in xrange(1, 3):
	t.mutate()
	t.pprint()
	print '\n\n'

print t.get_rand_sample(20)
print t.get_rand_sample(30)
print t.get_rand_sample(0)

t.print_img('tmp/')
