import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression
import argparse

parser = argparse.ArgumentParser(description='R2 values')
parser.add_argument('introns', type=str, help='file containing scored introns')
parser.add_argument('offset',type=int, help='offset')
parser.add_argument('slope',type=str, help='slope')
args = parser.parse_args()

f = open(args.introns)
exp = []
ime = []
while True:
    line = f.readline()
    if line == '': break

    data = line.split()
    exp.append(float(data[2][1:]))
    ime.append(float(data[-1]))


X = np.array(exp).reshape(-1,1)
y = ime
reg = LinearRegression().fit(X, y)
#print('Coefficients: ', reg.coef_)
r2 = str(reg.score(X,y))
txt = 'Offset: ' + str(args.offset) + ' Slope: ' + args.slope + ' R^2: '+r2 + '\n'
#print(txt)
f = open('Results','a')
f.write(txt)
f.close()
#y_pred = reg.predict((X))
#print(y_pred)

'''
plt.scatter(X,y, color='black')
plt.plot(X,y_pred, color='blue', linewidth=3)
plt.xticks(())
plt.yticks(())
plt.show()
'''
