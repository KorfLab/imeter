import math
import statistics
import sys


def kfreq(seqs, k, d, a):
	count = {}
	total = 0
	for s in seqs:
		for i in range(d, len(s) -k +1 -a):
			kmer = s[i:i+k]
			if kmer not in count: count[kmer] = 1
			count[kmer] += 1
			total += 1

	freq = {}
	for kmer in count: freq[kmer] = count[kmer] / total

	return freq



prox = []
dist = []
pdsplit = int(sys.argv[2])
k = int(sys.argv[3])
don = 5
acc = 10

with open(sys.argv[1]) as fp:
	for line in fp.readlines():
		f = line.split()
		beg = int(f[1])
		end = int(f[2])
		seq = f[-1]

		if end <= pdsplit: prox.append(seq)
		else:              dist.append(seq)

pfreq = kfreq(prox, k, don, acc)
dfreq = kfreq(dist, k, don, acc)

imeter = {}
for kmer in sorted(pfreq):
	imeter[kmer] = math.log2(pfreq[kmer] / dfreq[kmer])
#	print(f'{kmer}\t{pfreq[kmer]:.4f}\t{dfreq[kmer]:.4f}\t{math.log2(pfreq[kmer] / dfreq[kmer]):.4f}')

# score some introns
pscores = []
dscores = []
for seq in prox:
	score = 0
	for i in range(don, len(seq) -k +1 -acc):
		kmer = seq[i:i+k]
		score += imeter[kmer]
	pscores.append(score)

for seq in dist:
	score = 0
	for i in range(don, len(seq) -k +1 -acc):
		kmer = seq[i:i+k]
		score += imeter[kmer]
	dscores.append(score)

#print(pscores)
print(statistics.mean(pscores))
print(statistics.mean(dscores))
print(pscores)
