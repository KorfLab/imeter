
import argparse
import gzip
import itertools
import math
import sys

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

def linear(offset, pos, decay):
	if pos < offset: return 1.0
	w = 1 + (offset - pos) * decay
	if w < 0: return 0.0
	else:     return w

def geometric(offset, pos, decay):
	if pos < offset: return 1.0
	else:            return 1 / (1 + decay) ** (pos - offset)

def train(filename, model, offset, decay, k, a, d):

	# counts
	prox = kmers(k)
	dist = kmers(k)
	
	# weighting function
	if   model == 'linear':    func = linear
	elif model == 'geometric': func = geometric
	else:
		print(f'model {model} not implemented')
		sys.exit(1)

	if filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                        fp = open(filename)

	# counting
	for line in fp.readlines():
		f = line.split()
		beg = int(f[1])
		seq = f[-1]
		for i in range(d, len(seq) -k + 1 + a):
			kmer = seq[i:i+k]
			if kmer not in prox: continue
			x = beg + i
			pw = func(offset, x, decay)
			dw = 1 - pw
			prox[kmer] += pw
			dist[kmer] += dw

	# freqs
	pfreq = count2freq(prox)
	dfreq = count2freq(dist)
	imeter = {}
	for kmer in pfreq:
		imeter[kmer] = math.log2(pfreq[kmer] / dfreq[kmer])

	# done
	return imeter, pfreq, dfreq

# CLI

parser = argparse.ArgumentParser(description='IMEter version 3.0 trainer')
parser.add_argument('file', type=str, metavar='<file>',
	help='special intron file like at_ime_master.txt.gz')
parser.add_argument('--don', required=False, type=int, default=5,
	metavar='<int>', help='length of donor site [%(default)i]')
parser.add_argument('--acc', required=False, type=int, default=10,
	metavar='<int>', help='length of acceptor site [%(default)i]')
parser.add_argument('--model', required=False, type=str, default='linear',
	metavar='<name>', help='linear, geometric [%(default)s]')
parser.add_argument('--decay', required=False, type=float, default=0.01,
	metavar='<float>', help='decay constant [%(default).4f]')
parser.add_argument('--offset', required=False, type=int, default=100,
	metavar='<int>', help='decay offset [%(default)i]')
parser.add_argument('--k', required=False, type=int, default=5,
	metavar='<int>', help='kmer size [%(default)i]')
arg = parser.parse_args()

# Test
"""
for i in range(50):
	print(f'{i}\t{linear(10, i, 0.1):.3f}\t{geometric(10, i, 0.1):.3f}')
sys.exit(0)
"""

# Main

imeter, prox, dist = train(arg.file, arg.model, arg.offset, arg.decay,
	arg.k, arg.acc, arg.don)

# Output
print(f'# kmer, imeter, prox, dist ({arg.model} {arg.decay} {arg.offset})')
for k in imeter:
	print(k, imeter[k], prox[k], dist[k])
