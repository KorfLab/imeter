import re
import math
import statistics
import gzip
import itertools

def readintfile(infile):
	records = []
	if infile.endswith('.gz'): fp = gzip.open(infile, 'rt')
	else:                      fp = open(infile)

	for line in fp.readlines():
		f = line.split()
		pattern = '(\w+)\.1' #only look at the first entry, no splice variants
		match = re.search(pattern, f[0])
		if match:
			entry = {"Name":f[0], "beg":int(f[1]), "end":int(f[2]), "polarity":f[3],
			"Aerial":int(f[4]), "Carpel":int(f[5]), "DG Seed":int(f[6]),
			"LG Seed":int(f[7]), "Leaf":int(f[8]), "Pollen":int(f[9]),
			"Receptacle":int(f[10]), "RAM":int(f[11]), "Root":int(f[12]),
			"SAM":int(f[13]), "S12Flower":int(f[14]), "Gene":f[15]}

			records.append(entry)

	return(records)

def readmodelfile(infile):
	model = {}
	fp = None #lmao lexical scope

	if infile.endswith('.gz'): fp = gzip.open(infile, 'rt')
	else:                      fp = open(infile)

	for line in fp.readlines():
		line.rstrip()
		if line.startswith('$'): continue
		kmer, score, pfreq, kfreq = line.split()
		model[kmer] = float(score)
	k = len(list(model.keys())[0])\

	return(model, k)

def readfastafile(infile):
	genename = None
	fp = None
	seqs = []
	hold = []

	if infile.endswith('.gz'): fp = gzip.open(infile, 'rt')
	else:                      fp = open(infile)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(hold) > 0:
				seqs.append((name, ''.join(hold)))
				name = line[1:]
				hold = []
			else:
				name = line[1:]
		else:
			hold.append(line)

	fp.close()
	return(seqs)

def scoreintron(model, seq, k, d, a, cut):
	score = 0
	for i in range(d, len(seq) -k + 1 - a):
		kmer= seq[i:i+k]
		if kmer in model: score += model[kmer]
	return score

def generatekmers(k, init=0, alph='ACGT'):
	kmers = {}
	for tuple in itertools.product(alph, repeat=k):
		kmer= ''.join(tuple)
		kmers[kmer] = init

	return(kmers)

def decay(count, rate=0.002):
	#currently a linear decay after significance cutoff
	#rate of 0.002 will decay count to 0 after 1000 bp
	count = count - rate
	return(count)

def frequencies(counts):
	freqs = {}
	total = 0
	for kmer in counts: total += counts[kmer]
	for kmer in counts: freqs[kmer] = counts[kmer] / total
	return freqs

def train_imeter2(records, cut=400, k=5, d=5, a=10, sig=500):

	#initialize kmer counts
	prox = generatekmers(k)
	dist = generatekmers(k)
	count = 1

	#make counts
	for entry in records:
		seq = entry["Gene"]
		for i in range(d, len(seq) -k + 1 - a):
			kmer=seq[i:i+k]
			if kmer not in prox: continue
			if (entry["beg"] + i) >= sig:
				count = decay(count)
			if entry["beg"] <= cut: prox[kmer] += count
			else:                   dist[kmer] += count

	#convert to frequencies
	pfreqs = frequencies(prox)
	dfreqs= frequencies(dist)

	#make our training set
	imeter = {}

	for kmer in sorted(pfreqs):
		imeter[kmer] = math.log2(pfreqs[kmer]/dfreqs[kmer])

	return(imeter, pfreqs, dfreqs)
