import gzip
import math
import itertools
import sys
import argparse
#import kmer_weights

#using ArgParse
parser = argparse.ArgumentParser(description='IMEter trainer')
parser.add_argument('trainer', type=str, help='File used for getting training set')
args = parser.parse_args()

#Returns a dictionary of all kmers possible with a value of their count
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


D = 5 #length of splice donor site
K = 5 #kmer size
A = 10 #length of splice acceptor site
proximal_count = getAllKmers(K) #counts of each kmer - proximal
distal_count = getAllKmers(K) #counts of each kmer - distal

#get the counts
f = gzip.open(args.trainer, 'rt')
while True:
    line = f.readline()
    if line == '': break

    data = line.split()
    begin = data[1]
    strand = data[3]
    seq = data[-1]
    if strand != '+':
        continue

    for i in range(D,len(seq)-K-A+1):
        kmer = seq[i:i+K]
        if int(begin) < 400:
            if kmer in proximal_count:
                proximal_count[kmer] += 1
        else:
            if kmer in distal_count:
                distal_count[kmer] += 1
f.close()

#compute kmer frequencies
proximal_freq = getFrequency(proximal_count)
distal_freq = getFrequency(distal_count)

for kmer in proximal_freq:
    print(kmer, math.log2(proximal_freq[kmer]/distal_freq[kmer]))
