import gzip
import itertools
import math

prox_seqs = []
dist_seqs = []

t = 400
with gzip.open('../at_ime_master.txt.gz', 'rt') as fp:
	for line in fp.readlines():
		f = line.split()
		beg = int(f[1])
		seq = f[-1]
		if beg < 400: prox_seqs.append(seq)
		else:         dist_seqs.append(seq)

#print(len(prox_seqs), len(dist_seqs))

kcount = {}
ktotal = 0
P = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
total = 0
k = 5

for seq in prox_seqs:
	P['A'] += seq.count('A')
	P['C'] += seq.count('C')
	P['G'] += seq.count('G')
	P['T'] += seq.count('T')
	total += len(seq)

	for i in range(len(seq) -k +1):
		kmer = seq[i:i+k]
		if kmer not in kcount: kcount[kmer] = 0
		kcount[kmer] += 1
		ktotal += 1

for nt in P:
	P[nt] /= total

kfreq = {}
for kmer in kcount:
	kfreq[kmer] = kcount[kmer] / total

for t in itertools.product('ACGT', repeat=k):
	kmer = ''.join(t)
	p = 1
	for nt in kmer: p *= P[nt]

	print(kmer, kcount[kmer], kfreq[kmer], p,  p*total, math.log2(p*total/kcount[kmer]))
