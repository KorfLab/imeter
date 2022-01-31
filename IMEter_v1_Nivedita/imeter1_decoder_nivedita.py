import gzip
import korflib


D = 5 #length of splice donor site
K = 2 #kmer size
A = 10 #length of splice acceptor site

readfile = korflib.read_fasta('db_IME_Rose_WT_introns.fa')

print(readfile)
#f = gzip.open('mini.gz', 'rt')
'''
while True:
    line = f.readline()
    if line == '': break

    data = line.split()
    seq = data[-1]
    score = 0

    for i in range(D,len(seq)-K-A+1):
        kmer = seq[i:i+K]
        score += IMEter[kmer]
        print(kmer, score)
#f.close()
'''
