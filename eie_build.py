import glob
import sys
import grimoire.genome

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


## main ##
uid = 1
model = read_model('v1/imeter1.model') # we are very sorry for hard-coded paths

for chrom in grimoire.genome.Reader(fasta=sys.argv[1], gff=sys.argv[2]):
	xd = {} # expression data
	for f in chrom.ftable.features:
		if f.source == 'RNASeq_splice': xd[f'{f.beg}-{f.end}'] = f.score

	for gene in chrom.ftable.build_genes():
		if len(gene.transcripts()) == 0: continue # skip ncRNAs
		if gene.issues: continue # skip genes with obvious oddities
		for tx in gene.transcripts():
			if len(tx.introns) == 0: continue # skip intronless genes
			if len(tx.exons) != len(tx.cdss): continue # skip tx w/ utr exons
			for i, intron in enumerate(tx.introns):
				isig = f'{intron.beg}-{intron.end}'
				if isig not in xd: continue # skipping non-verified introns
				xs = int(xd[isig])
				if xs < 1000: continue # skip poorly expressed isoforms
				iseq = intron.seq_str()
				if len(iseq) > 500: continue # skip long introns			
				
				ime1 = score_intron(model, iseq, 5)
				if ime1 < 10: continue # skip low IME introns
				
				e5seq = tx.cdss[i].seq_str()
				e3seq = tx.cdss[i+1].seq_str()
				if len(e5seq) < 50: continue # need enough seq to align
				if len(e3seq) < 50: continue
				
				ib = intron.beg - gene.beg
				ie = intron.end - gene.beg
				if ib > 600: continue # prefer introns near TSS
				
				st = intron.strand
				rb = intron.beg - len(e5seq) # region begin
				re = intron.end + len(e3seq) # region end
				loc = f'{chrom.name}:{rb}-{re}'
				print(f'>eie-{uid} {loc} {st} {tx.id} {ib}-{ie} {xs} {ime1:.1f}')
				print(e5seq)
				print(iseq)
				print(e3seq)
				uid += 1


