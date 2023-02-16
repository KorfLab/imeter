# enumerate sequences based on their position weight matrix
# ie replace maybe a highish prob. of A or M
# decide using manhattan dist? 
"""
motif-finder - given a collection of sequences, find the motifs in common (like MEME)
motif-positioner - given a sequence and a motif, find all of the matching positions
motif-scorer - given a sequence and motif, provide an aggregate score or probability
mult-motif-scorer - as above, but with a collection of sequences and motifs
motif-comparer - given two motifs, measure how similar they are to each other
motif-displayer - it might be convenient to display digitized motifs as more traditional PWMs
"""

import numpy as np
from collections import Counter
import argparse
import math
import itertools

parser = argparse.ArgumentParser()
filename = parser.add_argument('filename') 

def parse_file(filename):
     with open(filename, 'r') as f:
        # read file into numpy array (rows: array of the sequence)
        lines = f.readlines()
        seq = [] 
        for line in lines:
            if not line.startswith('>') and line != "": 
                seq.append(np.array(list(line.strip())))
        seq = np.array(seq)
     return seq

def motif_positioner(sequences, motif):
    # finds all matching positions in format (index of sequence, index of motif in sequence)
    # remember 0 indexed!!
    k = len(motif)
    positions = []
    for ind, seq in enumerate(sequences):
        kmers = []
        for i in range(len(seq)-k+1):
            kmers.append(''.join(seq[i:i+k]))
        if motif in kmers:
            positions.append((ind, kmers.index(motif)))
    return positions

def count_kmers(seq, k, motif):
    # counts kmers
    kmers = []
    for i in seq:
        for j in range(len(i)-k + 1):
            kmers.append(''.join(i[j:j+k]))
    count_kmers = Counter(kmers)
    return count_kmers

def get_CG_content(sequences):
    # CG content for background probability
    unique, counts = np.unique(sequences, return_counts=True)
    output = {}
    tot = sum(counts)
    for u, c in zip(unique, counts):
        output[u] = c/tot
    return output

def prob_matrix(sequences):
    Bn = get_CG_content(sequences)
    prob = []
    tot = len(sequences.T[0])
    for seq in sequences.T: # iterate down columns
        freq = {'A':0, 'C':0, 'G':0, 'T':0}
        for i in seq:
            freq[i] += 1
        for key, value in freq.items(): # better way? broadcasting?
            freq[key] = value/tot # get prob.
        prob.append(freq) 
    return prob

def score_matrix(sequences):
    """
    score every letter:
    S(n) = log2[P(Cn)/Bn]
    Bn = background probability (use CG content or 0.25???)
    looks for prob of a certain letter in the seq being higher than that of the background prob 
    """
    scores = prob_matrix(sequences)
    Bn = get_CG_content(sequences)
    for i in scores:
        for key, value in i.items():
            if value > 0:
                i[key] = math.log2(value/Bn[key])
    return scores
    
def motif_comparer(motif_a, motif_b, seq):
    # compares a single column from each motif individually
    # and adds up the distance from all column comparisons.
    # Problems:
        # not sure if this is correct
        # calc. dist. from prob. of col. where motif matches?
        # also, not sure what to do if motif occurs more than once?
        # need pseudocounts so there's no zero probability??
        # what types of distance should we use (here I use euclidean)
    dist = 0
    matrix = prob_matrix(seq)
    for i in range(len(motif_a)):
        pos_a = motif_positioner(seq, motif_a)
        pos_b = motif_positioner(seq, motif_b)
        prob_motif_a = list(matrix[pos_a[0][1]].values())
        prob_motif_b = np.array(list(matrix[pos_b[0][1]].values()))
        dist += math.sqrt(sum(np.square((prob_motif_a - prob_motif_b))))
    return dist


def display():
    # display scores underneath like in fastq

    pass 


# QUESTION:
# What is a motif score????
# what is distance?
# comparing motifs:
     # need pseudocounts so there's no zero probability????
     # what to do if motif occurs more than once? 
# if using aligned sequences, what do we do with gaps?


def enumerate_seq(filename):
    # calc prob of each letter position by position
    # assign each position a letter
    seq = parse_file(seq)
    output_seq = ""
    for i in seq.T: # iterate down columns
        prob = Counter(i) # count freq of nuc. 
        # testing: using alphabet (ACGT)
        max_nuc = max(prob)
        # need to have rule for tiebreaker (currently by alphabetical order)
        if prob(max_nuc)/sum(prob.values()) < 0.5: # arbitrary?
            output_seq += 'N' # if prob. of max nucleotide < 50%, add N
        else:
            output_seq += max_nuc  # add nucleotide with max probability into the output sequence
    return output_seq