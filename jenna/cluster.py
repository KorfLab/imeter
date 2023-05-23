import pandas as pd
import seaborn as sns
from matplotlib import pyplot as plt
import sys

#Data set
df = pd.read_csv(sys.argv[1], sep=' ')
df = df.set_index('seq')


#default plot
sns.clustermap(df, metric ="euclidean", standard_scale=1, method="ward")


#show the graph
plt.show()
