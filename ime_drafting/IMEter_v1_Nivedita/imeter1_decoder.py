import gzip


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

def score(seq, D=5, K=2, A=10):
    score = 0
    for i in range(D,len(seq)-K-A+1):
        kmer = seq[i:i+K]
        score += IMEter[kmer]
        print(kmer, score)

'''
D = 5 #length of splice donor site
K = 2 #kmer size
A = 10 #length of splice acceptor site
'''
readfile = []
readfile = read_fasta('db_IME_Rose_WT_introns.fa')
for seq in readfile:
    print(score(seq[1]))


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
