import imelib
import sys

imeter, prox, dist = imelib.train_imeter1(sys.argv[1], k=3)

for k in imeter:
	print(k, imeter[k], prox[k], dist[k])
