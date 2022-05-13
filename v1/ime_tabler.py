import gzip

def read_model(filename):
	model = {}

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	for line in fp.readlines():
		if line.startswith('#'): continue
		kmer, imeter, prox, dist = line.split()
		model[kmer] = float(imeter)

	return model

def score_intron(model, seq, k, d=5, a=10):
	s = 0
	for i in range(d, len(seq) -k + 1 -a):
		kmer = seq[i:i+k]
		if kmer in model: s += model[kmer]
	return s


model = read_model('imeter1.model')
k = len(list(model.keys())[0])

with gzip.open('../at_ime_master.txt.gz', 'rt') as fp:
	for line in fp.readlines():
		f = line.split()
		name = f[0]
		beg = f[1]
		end = f[2]
		seq = f[-1]
		score = score_intron(model, seq, k)
		print(name, beg, end, score, sep=',')
