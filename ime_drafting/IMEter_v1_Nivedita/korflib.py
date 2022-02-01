import gzip
import random
import sys

def get_filepointer(filename):
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

def read_fasta(filename):

	name = None
	seqs = []

	fp = get_filepointer(filename)

	while True:
		line = fp.readline()
		if line == '': break
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

def revcomp(seq):
	comp = str.maketrans('ACGTRYMKWSBDHV', 'TGCAYRKMWSVHDB')
	anti = seq.translate(comp)[::-1]
	return anti

def randomseq(n):
	alph = 'ACGT'
	seq = []
	for i in range(n):
		seq.append(random.choice(alph))
	return ''.join(seq)


def parse_ab_blast(filename):
	
	fp = get_filepointer(filename)

	qid, sid = None, None
	qbeg, qend, sbeg, send = None, None, None, None
	score, pct, qstr, sstr = None, None, None, None
	data = False
	
	while True:
		line = fp.readline()
		if line == '':
			break
		elif line.startswith('Query='):
			f = line.split()
			qid, sid = f[1], None
			qbeg, qend, sbeg, send = None, None, None, None
			score, pct, qstr, sstr = None, None, None, None
			data = False
		elif line.startswith('Parameters:'):
			if data:
				yield qid, sid, qbeg, qend, sbeg, send, score, pct, qstr, sstr
				qbeg, qend, sbeg, send = None, None, None, None
				score, pct, qstr, sstr = None, None, None, None
				data = False
		elif line.startswith('>'):
			if data:
				yield qid, sid, qbeg, qend, sbeg, send, score, pct, qstr, sstr
				qbeg, qend, sbeg, send = None, None, None, None
				score, pct, qstr, sstr = None, None, None, None
				data = False
			f = line.split()
			sid = f[0][1:]
		elif line.startswith(' Score ='):
			if data:
				yield qid, sid, qbeg, qend, sbeg, send, score, pct, qstr, sstr
				qbeg, qend, sbeg, send = None, None, None, None
				score, pct, qstr, sstr = None, None, None, None
				data = False
			data = True
			f = line.split()
			score = int(f[2])
			line2 = fp.readline()
			f = line2.split()
			ids = f[2]
			match, total = ids.split('/')
			pct = int(match)/int(total)
			qstr, sstr = '', ''
		elif line.startswith('Query:'):
			f = line.split()
			qstr += f[2]
			if qbeg == None: qbeg = f[1]
			qend = f[3]
		elif line.startswith('Sbjct:'):
			f = line.split()
			sstr += f[2]
			if sbeg == None: sbeg = f[1]
			send = f[3]

	fp.close()
	

def parse_ab_blast_mformat(filename):
	
	fp = get_filepointer(filename)

	qid, sid = None, None
	qbeg, qend, sbeg, send = None, None, None, None
	score, pct, qstr, sstr = None, None, None, None
	data = False
	
	while True:
		line = fp.readline()
		if line == '': break
		if line.startswith('#'): continue
		qid, sid, E, N, Sp, score, alen, nid, npos, nsim, pct, pct2,\
			qgn, qqlen, sgn, sqlen, qf, qbeg, qend, sf, sbeg, send\
			= line.split()
		yield qid, sid, int(qbeg), int(qend), int(sbeg), int(send),\
			float(score), float(pct)
	
	fp.close()
