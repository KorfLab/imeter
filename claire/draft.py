import numpy as np
from collections import Counter
import argparse


def p_value(sequences):
    # p-value of proximal vs distal introns
    # 1. Get prob that kmer is proximal?
    # 2. divide by 
    prox, dist = parse_file(filename)
    pkmers = get_kmers(prox, 4)
    dkmers = get_kmers(dist, 4)
    prob_prox = get_CG_content(prox)
    

    
def get_kmers(sequences, k):
    x = []
    for row in sequences:
        for col in range(len(row)-k+1):
            x.append(row[i:i+k])
    kmers = Counter(x)
    return kmers
    
    
def get_CG_content(sequences):
    # could also do for dinuc.
    n, m = np.shape(sequences)
    tot = n * m
    freq = {'A':0,'T':0,'C':0, 'G':0}
    for col in sequences.T:
        for row in col:
            freq[row] += 1
    for key, value in freq:
        freq[key] = value/tot
    return freq
        
        
    
def parse_file(filename):
    dist_sequences = []
    prox_sequences = []
    with gzip.open(filename, "rt") as f:
        lines = f.readlines()
        proximal = {}
        distal = {}
        for line in lines:
            trans_id, start, end, strand, exp, seq = line.strip().split()
            if start < 400:
                prox_sequences.append(np.array(seq.split()))
            else:
                dist_sequences.append(np.array(seq.split()))

    prox_sequences = np.array(sequences)
    dist_sequences = np.array(sequences)
    return prox_sequences, dist_sequences

# expected kmers per seq (1 per or many on one)

# getting p-values might be hard b/c numbers are large
# what's the most interesting (absolute vs relative magnitude)
# 10 vs 0 or 100:10
# poisson

# get prob. models 
# optimize search (most interesting kmers, pvalues?)