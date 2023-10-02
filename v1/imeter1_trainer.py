
import argparse
import gzip
import itertools
import math


def kmers(k, init=0, alph='ACGT'):
	kmers = {}
	for tup in itertools.product(alph, repeat=k):
		kmer = ''.join(tup)
		kmers[kmer] = init
	return kmers

def count2freq(count):
	freq = {}
	total = 0
	for k in count: total += count[k]
	for k in count: freq[k] = count[k] / total
	return freq

def train_imeter1(filename, k=5, d=5, a=10, t=400):

	# counts
	prox = kmers(k, init=1)
	dist = kmers(k, init=1)

	if filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                        fp = open(filename)

	for line in fp.readlines():
		f = line.split()
		beg = int(f[1])
		seq = f[-1]
		for i in range(d, len(seq) -k + 1 - a):
			kmer = seq[i:i+k]
			if kmer not in prox: continue
			if beg <= t: prox[kmer] += 1
			else:        dist[kmer] += 1

	# freqs
	pfreq = count2freq(prox)
	dfreq = count2freq(dist)
	imeter = {}
	for kmer in pfreq:
		imeter[kmer] = math.log2(pfreq[kmer] / dfreq[kmer])

	# done
	return imeter, pfreq, dfreq

# CLI

parser = argparse.ArgumentParser(description='IMEter version 1.0 trainer')
parser.add_argument('file', type=str, metavar='<file>',
	help='special intron file like at_ime_master.txt.gz')
parser.add_argument('--don', required=False, type=int, default=5,
	metavar='<int>', help='length of donor site [%(default)i]')
parser.add_argument('--acc', required=False, type=int, default=10,
	metavar='<int>', help='length of acceptor site [%(default)i]')
parser.add_argument('--cut', required=False, type=int, default=400,
	metavar='<int>', help='proximal-distal cutoff [%(default)i]')
parser.add_argument('--k', required=False, type=int, default=5,
	metavar='<int>', help='kmer size [%(default)i]')
arg = parser.parse_args()

# Main

imeter, prox, dist = train_imeter1(arg.file, t=arg.cut, k=arg.k,
	a=arg.acc, d=arg.don)

# Output
print('# kmer, imeter, prox, dist')
for k in imeter:
	print(k, imeter[k], prox[k], dist[k])
