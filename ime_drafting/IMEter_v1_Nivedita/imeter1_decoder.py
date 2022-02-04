import gzip
import sys

def get_filepointer(filename):
	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	elif filename == '-':          fp = sys.stdin
	else:                          fp = open(filename)
	return fp

def read_fasta(filename):

	name = None
	seqs = []

	fp = get_filepointer(filename)

	while True:
		line = fp.readline()
		if line == '': break
		line = line.rstrip()
		if line.startswith('>'):
			if len(seqs) > 0:
				seq = ''.join(seqs)
				yield(name, seq)
				name = line[1:]
				seqs = []
			else:
				name = line[1:]
		else:
			seqs.append(line)
	yield(name, ''.join(seqs))
	fp.close()

def imeter(seq, K, model, D=5, A=10):
    s = 0
    for i in range(D,len(seq)-K-A+1):
        kmer = seq[i:i+K]
        s += model[kmer]
    return s

model = {}
K = None
fp = open(sys.argv[1])
for line in fp.readlines():
    kmer, score = line.split()
    K = len(kmer)
    model[kmer] = float(score)



for name, seq in read_fasta('db_IME_Rose_WT_introns.fa'):
    s = imeter(seq, K, model)
    print(name, s)


#f = gzip.open('mini.gz', 'rt')
'''
def score(D,):
    for i in range(D,len(seq)-K-A+1):
        kmer = seq[i:i+K]
        score += IMEter[kmer]
        print(kmer, score)


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
