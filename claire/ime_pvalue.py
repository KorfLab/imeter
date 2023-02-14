import numpy as np
import argparse
import gzip

parser = argparse.ArgumentParser()
parser.add_argument("filename", help="enter read file", type=str)
parser.add_argument("k", help="kmer length", type=int)
args = parser.parse_args()

def parse_file(filename):
    # read through file and add sequences to numpy array
    dist_sequences = []
    prox_sequences = []
    with gzip.open(filename, "rt") as f:
        lines = f.readlines()
        training_dataset = lines[0:int(len(lines)/4)]
        testing_dataset = lines[int(len(lines)/4):]
    return training_dataset, testing_dataset
    
def get_prox_dist(sequences):
    prox_sequences = []
    dist_sequences = []
    for line in sequences:
        line = line.rstrip().split()
        start = int(line[1])
        end = int(line[2])
        seq = line[15]
        if start < 400:
            prox_sequences.append(np.array([*seq]))
        else:
            dist_sequences.append(np.array([*seq]))  
    prox_sequences = np.array(prox_sequences)
    dist_sequences = np.array(dist_sequences)
    return prox_sequences, dist_sequences
    
def get_kmers(sequences, k):
    x = []
    for row in sequences:
        for col in range(len(row)-k+1):
            x.append(list(row[col:col+k]))
    kmers = []
    for i in x:
        kmers.append("".join(i))
    return kmers
    
    
def get_CG_content(sequences):
    # could also do for dinuc.
    freq = {'A':0,'T':0,'C':0, 'G':0}
    for col in sequences.T:
        for row in col:
            freq[row] += 1
    tot = sum(freq.values())
    for key, value in freq.items():
        freq[key] = value/tot
    return freq
    
def p_value(sequences):
    # p-value of proximal vs distal introns
    # z = (y - mu)/sigma
    # mu = avg = average %T in proximal?
    # sigma = std dev of %T in proximal
    # maybe using R, test multiple regression with categorical variables (including which % of nuc predicts prox / distal the best
    # x's: T% (dependent on other nuc.? b/c they add up to 1)
    # y's: categorical
    prox, dist = parse_file(filename)
    pkmers = get_kmers(prox, 4)
    dkmers = get_kmers(dist, 4)
    prob_prox = get_CG_content(prox)
    expected = 0

training_dataset, testing_dataset = parse_file(args.filename) 
prox_sequences, dist_sequences = get_prox_dist(training_dataset)
prox_CG = get_CG_content(prox_sequences)
dist_CG = get_CG_content(dist_sequences)
print(prox_CG, dist_CG)



# expected kmers per seq (1 per or many on one)

# getting p-values might be hard b/c numbers are large
# what's the most interesting (absolute vs relative magnitude)
# 10 vs 0 or 100:10
# poisson

# get prob. models 
# optimize search (most interesting kmers, pvalues?)

### Logistic regression for classification
"""
#Preparing the model
from sklearn.model_selection import train_test_split
X_train, X_test, y_train, y_test = train_test_split(X,y, test_size = 0.20, random_state = 99)
from sklearn.linear_model import LogisticRegression
model = LogisticRegression()
model.fit(X_train, y_train)
model.score(X_test,y_test)
"""