import math
import statistics
import sys
import gzip
import itertools as iter
import csv

def generatefilters(): #generates a list of legal kmers. in the future, user specified exceptions could be handled
    filter = {}
    perms = list(iter.product('ACTG', repeat=5))
    for i in range(0, len(perms)):
        perms[i] = ''.join(perms[i])
        filter[perms[i]] = 1
    return(filter)

def readinfile(filename): #Reads the entire file into memory.
    records = []
    with gzip.open(filename, 'rt') as fp:
        for line in fp.readlines():
            f = line.split()
            entry = {'Name':f[0], 'beg':int(f[1]), 'end':int(f[2]), 'polarity':f[3], \
            'Aerial':int(f[4]), 'Carpel':int(f[5]), 'DG Seed':int(f[6]), 'LG Seed':int(f[7]), 'Leaf':int(f[8]), \
            'Pollen':int(f[9]), 'Receptacle':int(f[10]), 'RAM':int(f[11]), 'Root':int(f[12]), 'SAM':int(f[13]), \
            'S12Flower':int(f[14]), 'Gene':f[15]} #this is an unholy abomination. abandon hope, all ye who enter here.
            records.append(entry)
    fp.close()
    return(records)

def cutoffsplit(records, cutoff): #splits sequences into proximal and distal based on start
    prox = []
    dist = []
    for i in range(0, len(records)):
        if records[i]['beg'] <= cutoff:
            prox.append(records[i]['Name'], records[i]['Gene'])
        else:
            dist.append(records[i]['Name'], records[i]['Gene'])
    return(prox, dist) #returns sequences split into two

def kmercount(seqs, k): #determines the total count of kmers in the sequences
    don = 5 #donor seqeuence
    acc = 10 #acceptor sequence. we're...hardcoding these, I guess. Not a fan. Will these ever need to change?
    total = 0
    decay = 1
    count = {}
    filter = generatefilters() #a dictionary of all the possible legal kmers. reformatting to dictionary drastically cuts down processing time

    for s in seqs:
        for i in range(don, len(s[1])-k+1-acc):
            kmer = s[1][i:i+k]
            if kmer in filter: #this is a CPU sink, really big O
                if kmer not in count:
                    count[kmer] = 0
                count[kmer] += decay #to factor in geom decay, to account for intron significance dropping off as it lengthens
                total += decay #pretty sure we want to make the total with the decay as well
                #some sort of decay function goes here...

    return(count, total) #returns a tuple, a dictionary keyed with unique kmers and their counts, and the total number)

def kmerfreqs(counts): #determines frequencies of kmers across an introns
    #print(counts)
    freqs= {}
    for kmer in counts[0]:
        freqs[kmer] = counts[0][kmer] / counts[1]
    return(freqs)

def training(proxfreqs, distalfreqs, xfold): #calculates the log odd probability for our kmerfreqs
    trained = {}
    for kmer in sorted(proxfreqs):
        trained[kmer] = math.log2(proxfreqs[kmer] / distalfreqs[kmer])
        #write out training set: kmer, proxfreq, distalfreq, score
    with open("trainedimeter.txt", "w", newline='') as outcsv:
        trwriter = csv.writer(outcsv, delimiter=',', quotechar = '"', quoting = csv.QUOTE_NONE, escapechar='|')
        trwriter.writerow(['kmer','proximal frequency','distal frequency','log-odds score'])
        for kmer in sorted(proxfreqs):
            trwriter.writerow([kmer,proxfreqs[kmer],distalfreqs[kmer],trained[kmer]])
    outcsv.close()

    return(trained)

def scoring(prox, distal, records, trained, k): #calculates the score of a query. Possible support for dynamic queries
    proxscores = []
    distalscores = []
    don = 5
    acc = 10

    for seq in prox:
        score = 0
        for i in range(don, len(seq) -k +1 -acc):
            kmer = seq[i:i+k]
            if kmer in trained:
                score += trained[kmer]
        proxscores.append((seq, score))

    for seq in distal:
        score = 0
        for i in range(don, len(seq) -k +1 -acc):
            kmer = seq[i:i+k]
            if kmer in trained:
                score += trained[kmer]
        distalscores.append((seq, score))
    return(proxscores, distalscores)


def printscores(records, outfile):
    if outfile is None:
        print('f')
