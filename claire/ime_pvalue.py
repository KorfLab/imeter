import numpy as np
import argparse
import gzip
import statistics

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
    # get avg / std dev of nuc. % dist.
    # could also do for dinuc.
    freq = {'A':0,'T':0,'C':0, 'G':0}
    freqs = {'A':[],'T':[],'C':[], 'G':[]}
    stdDev = {}
    avg = {}
    for ind, row in enumerate(sequences): 
        for col in row:
            freq[col] += 1
        for key in freqs.keys():
            freqs[key].append(freq[key])
        freq = {'A':0,'T':0,'C':0, 'G':0}
        tot = len(row)
        for key, value in freqs.items(): # get probability
            freqs[key] = value[:ind] + [value[ind:ind+1][0]/tot] + value[ind+1:] 
    for key, value in freqs.items():
        avg[key] = sum(value)/len(value)
        stdDev[key] = statistics.stdev(value)   
    return freqs, avg, stdDev
    
def predicting_prox_dist(prox_seq, dist_seq):
    # naive bayes:
        # p(C|x) = p(C)p(x|C) / p(x)
        # x = feature / trait, C = category
    # Gaussian naive bayes
    # prob density 
    # P(x = trait | category) = 1/ sqrt(2*pi*var) * exp(-(trait - mean)^2 / (2*var))
    pfreqs, pavg, pstdDev = get_CG_content(prox_sequences)
    dfreqs, davg, dstdDev = get_CG_content(dist_sequences)
    prob_C = len(prox_seq) / ((len(prox_seq) + len(dist_seq)
    
    
def p_value(prox, distsequences):
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
    mean_prox, stdev_prox = get_CG_content(prox)
    mean_dist, stdev_dist = get_CG_content(dist)
    

training_dataset, testing_dataset = parse_file(args.filename) 
prox_sequences, dist_sequences = get_prox_dist(training_dataset)
prox_CG, prox_avg, prox_stdDev = get_CG_content(prox_sequences)
dist_CG, dist_avg, dist_stdDev = get_CG_content(dist_sequences)
print("Proximal average: ", prox_avg)
print("Proximal std dev: ", prox_stdDev)
print("Distal average: ", dist_avg)
print("Distal std dev: ", dist_stdDev)



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