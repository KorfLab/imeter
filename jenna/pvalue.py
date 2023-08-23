import gzip

with gzip.open('../at_ime_master.txt.gz', 'rt') as fp:
	for line in fp.readlines():
		f   = line.split()
		beg = f[1]
		seq = f[-1]
		print(beg, seq)
