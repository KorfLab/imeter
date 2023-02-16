import numpy as np
import regex as re
import gzip
import argparse


parser = argparse.ArgumentParser()
filename = parser.add_argument('filename') 
args = parser.parse_args()


def read_file(filename):
     with gzip.open(filename, 'r') as f:
        # read file into numpy array (rows: array of the sequence)
        lines = f.readlines()
        seq = [] 
        for line in lines:
            if not line.startswith('>') and line != "": 
                seq.append(np.array(list(line.strip())))
        seq = np.array(seq)
     return seq
     
def read_jaspar(filename):
    freq = []
    prob = []
    with open(filename) as f:
        lines = f.readlines()
        for line in lines:
            line = line.rstrip()
            if not line.startswith('>'):
                freq.append(np.array(re.findall("\d+", line))) # get digits (rows = A,C,G,T)
    freq = np.array(freq)
    x = []
    for col in freq.T:
        for row in col:
            x.append(int(row))
        x = [i/sum(x) for i in x]
        prob.append(x)
        x = []
    return prob


def reverse_complement(motif):
	nuc = {'A':'T', 'C':'G', 'T':'A', 'G':'C'}
	rc = ''
	#for key in 
	#return out
	

	

def motif_comparer(motif1, motif2):
	# compare using manhattan distance (|x2-x1|)
	# assuming motifs have the same len
	# input: motifs with prob. dist. (dict) at each position
	score = []
	tscore = 0
	for i, j in enumerate(motif1, motif2):
		for key in motif1.keys():
			tscore = math.abs(motif1[key] - motif2[key])
		score.append(tscore)
		tscore = 0
	return score

freq = read_jaspar(args.filename)
print(freq)
#motifA =  [{'A':0.5, 'T':0.5, 'C':0, 'G':0}, {'A':0, 'T':0.5, 'C':0.25, 'G':0.25}, {'A':0, 'T':0, 'C':1, 'G':0}, {'A':0.25, 'T':0, 'C':0, 'G':0.75}]
#motifB =  [{'A':0.5, 'T':0.5, 'C':0, 'G':0}, {'A':0, 'T':0.5, 'C':0.25, 'G':0.25}, {'A':0, 'T':0, 'C':0, 'G':1}, {'A':0, 'T':0.3, 'C':0, 'G':0.7}]
motifA = [

	[.5,.5,0,0],
	[0,.3,.7,0],
	[.2,.2,.2,.4],
	[.9,.1,0,0]
]

# output file format header: % model_type name length
# sample prob: 



