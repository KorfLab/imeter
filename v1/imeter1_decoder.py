
import argparse
import math

def read_fasta(filename):
	name = None
	seqs = []

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	for line in fp.readlines():
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

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

# CLI

parser = argparse.ArgumentParser(description='IMEter version 1.0 decoder')
parser.add_argument('model', type=str, metavar='<file>',
	help='IMEter model file')
parser.add_argument('fasta', type=str, metavar='<file>',
	help='fasta file')
arg = parser.parse_args()

# Main

model = read_model(arg.model)
k = len(list(model.keys())[0])

for desc, seq in read_fasta(arg.fasta):
	name = desc.split(maxsplit=1)[0]
	print(name, score_intron(model, seq, k))

