#!/usr/bin/env python3

import sys
import os

from grimoire.genome import Reader

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

# Main (python3 intron_report.py TAIR10_chr_all.fas TAIR10_GFF3_genes.gff ../v1/imeter1.model > table)

model = read_model(sys.argv[3])
ime_k = len(list(model.keys())[0])

genome = Reader(fasta=sys.argv[1], gff=sys.argv[2])
done = {} # to remove redundancy
for chrom in genome:
	for gene in chrom.ftable.build_genes():
		if len(gene.transcripts()) == 0: continue # skip ncRNAs
		for tx in gene.transcripts():
			if len(tx.introns) == 0: continue # skip intronless genes
			for intron in tx.introns:
				iseq = intron.seq_str()
				score = score_intron(model, iseq, ime_k)
				sig = (chrom.name, intron.beg, intron.end)
				if sig in done: continue # redundant
				done[sig] = True
				print(chrom.name, intron.beg, intron.end, intron.strand, score,
					tx.id)

