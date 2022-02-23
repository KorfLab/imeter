import matplotlib.pyplot as plt
import numpy as np
from sklearn import datasets, linear_model
from sklearn.metrics import mean_squared_error, r2_score
from sklearn.linear_model import LinearRegression


X = np.array([[12, 6], [10, 5], [4, 2], [0, 0]])
y = np.dot(X, np.array([1, 2])) + 3
#ime = [12, 10, 4, 0]
#exp = [6, 5, 2, 0]
regr = linear_model.LinearRegression()
reg = LinearRegression().fit(X, y)
print(reg.score(X, y))
