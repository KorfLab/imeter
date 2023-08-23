import math
import sys
import statistics
import gzip

def shape(vals):
	s = sum(vals)
	if s == 0:
		return None
	p = []
	for val in vals:
		p.append(val/s)
	return entropy(p)

def entropy(ps):
	h = 0
	for p in ps:
		if p != 0:
			h -= p*math.log2(p)
	return h

def read_model(filename):
	model = {}

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	for line in fp.readlines():
		if line.startswith('#'): continue
		kmer, imeter = line.split()
		model[kmer] = float(imeter)

	return model

def score_intron(model, seq, k, d=5, a=10):
	s = 0
	for i in range(d, len(seq) -k + 1 -a):
		kmer = seq[i:i+k]
		if kmer in model: s += model[kmer]
	return s
	
	
def categorize(bounds, val):
	for i,b in enumerate(bounds):
		if val < b: return i
	return len(bounds) 
	
ime_model = read_model("dm_ime_model")

ibs   = [200, 400, 800, 1600]
ds    = [100, 200, 400, 800]
imes  = [-10, 0, 10, 20, 40]
prox  = []
dist  = []
c_scores = {}
scores = {}

fp = open("../datacore/proteomes/fly_v_worm")
for line in fp.readlines():
	f = line.split()
	name   = f[0]
	pct_id = f[2]
	if name not in scores: scores[name] = [pct_id]
	

fp = open("../datacore/proteomes/dm_names.txt")
for line in fp.readlines():
	f = line.split()
	pid = f[1]
	tid = f[2]
	if pid in scores: scores[pid].append(tid) 

tscores = scores.values()

for c, tid in tscores:
	if tid not in c_scores: c_scores[tid] = c

with gzip.open(sys.argv[1], "rt") as fp:
	for line in fp.readlines():
		line = line.rstrip()
		f = line.split()
		x = float(f[4]) 
		if x == 0: continue
		x = math.log2(x)
		ib = categorize(ibs,int(f[1]))
		ie = int(f[2])
		d = categorize(ds, (int(f[2]) - int(f[1])))
		ime = categorize(imes, score_intron(ime_model, f[-1], 4))
		if f[0] in c_scores: c = c_scores[f[0]]
		else: c = 0
		print(f[0], x, ib, d, ime, c)
		
		
'''
categories:
	ib: (200, 400, 800, 1600)
	d: (100, 200, 400, 800)
	ime: (-10, 0, 10, 20, 40)
	
	dm: 
		median x: ~2000
		top 88.66%: 50939 cutoff for 10000 x
		
'''
