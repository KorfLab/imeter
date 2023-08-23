import sys
import random

def read_model(filename):
	model = {}

	fp = None
	if   filename.endswith('.gz'): fp = gzip.open(filename, 'rt')
	else:                          fp = open(filename)

	for line in fp.readlines():
		if line.startswith('#'): continue
		kmer, imeter, prox, dist = line.split()
		model[kmer] = float(imeter)

	return model

def score_intron(model, seq, k, d=5, a=10):
	s = 0
	for i in range(d, len(seq) -k + 1 -a):
		kmer = seq[i:i+k]
		if kmer in model: s += model[kmer]
	return s

def fitness(model, seq):
	gc = (seq.count('G') + seq.count('C')) / len(seq)
	ime = score_intron(model, seq, 5)
	score = ime * (1 - abs(0.35 - gc)) ** 3
	return score

def create_individual(size, model):
	seq = ''
	for i in range(size): seq += random.choice('ACGT')
	return seq, fitness(model, seq)

def mate(p1, p2, mf, model):

	# recombination
	k = 5
	seq = ''
	while True:
		if random.random() < 0.5: parent = p1
		else: parent = p2
		pos = random.randint(0, len(parent[0]) -k + 1)
		seq += parent[0][pos:pos+k]
		if len(seq) > len(parent[0]):
			seq = seq[:len(parent[0])]
			break
	
	# mutation
	letters = list(seq)
	for i in range(len(letters)):
		if random.random() < mf:
			letters[i] = random.choice('ACGT')
	seq = ''.join(letters)
	
	return seq, fitness(model, seq)
	

# CLI
POPSIZE = 2000
SEQLEN = 300
GENERATIONS = 100
MUTATION = 0.01

# init
model = read_model(sys.argv[1])
#random.seed(1)

# create an initial population of individuals with random genotypes
population = []
for i in range(POPSIZE):
	population.append(create_individual(SEQLEN, model))

# breed population for n generations
for i in range(GENERATIONS):
	population.sort(key = lambda x: x[1], reverse=True)
	print(i, population[0][1], population[-1][1], population[0][0][10:40])
	for j in range(POPSIZE//2, POPSIZE):
		p1 = random.choice(population[:POPSIZE//2])
		p2 = random.choice(population[:POPSIZE//2])
		population[j] = mate(p1, p2, MUTATION, model)

best = population[0][0]
g = best.count('G')
c = best.count('C')
print((g+c)/SEQLEN)

print(population[0][0])