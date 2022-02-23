import gzip
import math
import itertools
import sys
import argparse

#using ArgParse
parser = argparse.ArgumentParser(description='IMEter trainer')
parser.add_argument('trainer', type=str, help='File used for getting training set')
args = parser.parse_args()

#Returns a dictionary of all kmers possible with a value of their count set to 0
def getAllKmers(size):
    allKmers = {}
    for kmer in itertools.product('ACGT', repeat=size):
        kmer = ''.join(kmer)
        allKmers[kmer] = 0
    return allKmers

#Returns a dictionary that keeps track of the kmer and its corresponding frequency
def getFrequency(kmers):
    total = 0
    for k in kmers:
        total += kmers[k]
    frequency = {}
    for kmer in kmers:
        frequency[kmer] = kmers[kmer] / total
    return frequency

#Assigns different weights based on where the kmer is???
def linear(num, start, slope):
    if num < start:
        p_weight = 1.0
        d_weight = 0.0
    elif (num - start) * slope >= 1:
        p_weight = 0.0
        d_weight = 1.0
    else:
        d_weight = (num - start) * slope
        p_weight = 1 - d_weight
    return p_weight, d_weight

D = 5 #length of splice donor site
K = 5 #kmer size
A = 10 #length of splice acceptor site
proximal_count = getAllKmers(K) #counts of each kmer - proximal
distal_count = getAllKmers(K) #counts of each kmer - distal

# Do we need proximal/distal if we are using weights?

#get the counts
f = gzip.open(args.trainer, 'rt')
while True:
    line = f.readline()
    if line == '': break

    data = line.split()
    #begin = data[1]
    strand = data[3]
    seq = data[-1]
    if strand != '+':
        continue

    for i in range(D,len(seq)-K-A+1):
        kmer = seq[i:i+K]
        if kmer in proximal_count:
            # counts proximal value according to weight
            proximal_count[kmer] += linear(i, 400, 1/400)[0]
        if kmer in distal_count:
            # counts distal value according to weight
            distal_count[kmer] += linear(i, 400, 1/400)[1]
        #print(kmer, proximal_count[kmer], distal_count[kmer])
f.close()

#compute kmer frequencies
proximal_freq = getFrequency(proximal_count)
distal_freq = getFrequency(distal_count)

for kmer in proximal_freq:
    print(kmer, math.log2(proximal_freq[kmer]/distal_freq[kmer]))
