import numpy as np
from sklearn.metrics import r2_score
from sklearn.linear_model import LinearRegression


X = np.array([[12, 6], [10, 5], [4, 2], [0, 0]])
y = np.dot(X, np.array([1, 2])) + 3
#ime = [12, 10, 4, 0]
#exp = [6, 5, 2, 0]
reg = LinearRegression().fit(X, y)
print(reg.score(X, y))
