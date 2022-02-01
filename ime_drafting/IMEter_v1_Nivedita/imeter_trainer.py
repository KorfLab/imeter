import gzip
import math
import itertools
import korflib

#Returns a dictionary of all kmers possible with a value of their count
def getAllKmers(size):
    allKmers = {}
    for kmer in itertools.product('ACGT', repeat=size):
        kmer = ''.join(kmer)
        allKmers[kmer] = 0
    return allKmers

#Returns a dictionary that keeps track of the kmer and its corresponding frequency
def getFrequency(kmerList, kmerCount):
    frequency = {}
    for kmer in kmerList:
        frequency[kmer] = kmerList[kmer] / kmerCount
    return frequency

#Returns a dictionary that contains all the kmers and their IMEter scores
def IMEter(proxFrequency,distFrequency):
    IMEter = {}
    for kmer in proxFrequency:
        IMEter[kmer] = math.log2(proxFrequency[kmer]/distFrequency[kmer])
    return IMEter

D = 5 #length of splice donor site
K = 2 #kmer size
A = 10 #length of splice acceptor site
totalprox = 0 #total # of proximal kmers
totaldist = 0 #total # of distal kmers
proximal_count = getAllKmers(K) #counts of each kmer - proximal
distal_count = getAllKmers(K) #counts of each kmer - distal

# get the counts

f = gzip.open('mini.gz', 'rt')
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
                totalprox += 1
        else:
            if kmer in distal_count:
                distal_count[kmer] += 1
                totaldist += 1
f.close()

#compute kmer frequencies
proximal_freq = getFrequency(proximal_count,totalprox)
distal_freq = getFrequency(distal_count,totaldist)

#calculate IMEter
IMEter = IMEter(proximal_freq, distal_freq)
