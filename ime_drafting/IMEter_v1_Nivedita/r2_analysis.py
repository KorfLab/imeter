import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression
import argparse

parser = argparse.ArgumentParser(description='R2 values')
parser.add_argument('introns', type=str, help='file containing scored introns')
args = parser.parse_args()

f = open(args.introns)
exp = []
ime = []
while True:
    line = f.readline()
    if line == '': break

    data = line.split()
    #print(data)
    exp.append(float(data[2][1:]))
    ime.append(float(data[-1]))


X = np.array(exp).reshape(-1,1)
y = ime
#a = np.array(ime).reshape(-1,1)
regr = linear_model.LinearRegression()
#print(X,y)
reg = LinearRegression().fit(X, y)
print(reg.score(X,y))
#print(reg.predict(a))
