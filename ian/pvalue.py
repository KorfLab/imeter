import gzip

prox = []
dist = []

t = 400
with gzip.open('../at_ime_master.txt.gz', 'rt') as fp:
	for line in fp.readlines():
		f = line.split()
		beg = int(f[1])
		seq = f[-1]
		if beg < 400: prox.append(seq)
		else:         dist.append(seq)

print(len(prox), len(dist))
