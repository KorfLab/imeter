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

def categorize(bounds, val):
	for i,b in enumerate(bounds):
		if val < b: return i
	return len(bounds) 

ibs  = [200, 400, 800, 1600]
ds   = [100, 200, 400, 800]
prox = []
dist = []

c_scores = {}
f_scores = []
scores = {}

fp = open("../datacore/proteomes/worm_v_fly")
for line in fp.readlines():
	f = line.split()
	name   = f[0]
	pct_id = f[2]
	if name not in scores: scores[name] = [pct_id]

fp = open("../datacore/proteomes/ce_names.txt")
for line in fp.readlines():
	f = line.split()
	pid = f[0]
	tid = f[2]
	if pid in scores: scores[pid].append(tid) 

t_scores = scores.values()

for s in t_scores:
	if len(s) == 2: f_scores.append(s)

for c, tid in f_scores:
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
		name = f[0][11:]
		if name in c_scores: c = c_scores[name]
		else: c = 0
		print(name, x, ib, d, c)
		
		
'''
categories:
	ib: (200, 400, 800, 1600)
	d: (100, 200, 400, 800)
	ime: (-10, 0, 10, 20, 40)		
'''
