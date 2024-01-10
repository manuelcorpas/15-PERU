import pandas as pd
from sklearn.manifold import MDS
from scipy.cluster.hierarchy import linkage, dendrogram
import matplotlib.pyplot as plt

# Load D-statistics data
file_path = 'ANALYSIS/15-DI-STAT/d_statistics_results.txt'
df = pd.read_csv(file_path, delim_whitespace=True)

# Creating a distance matrix from D-statistics
distance_matrix = df.pivot_table(index='pop2', columns='pop3', values='est', aggfunc=lambda x: abs(x.mean()))

# Make the matrix square
all_pops = set(distance_matrix.index).union(distance_matrix.columns)
distance_matrix = distance_matrix.reindex(index=all_pops, columns=all_pops)

# Explicitly symmetrize the matrix
symmetrized_matrix = (distance_matrix + distance_matrix.T) / 2

# Fill any remaining NaN values
symmetrized_matrix.fillna(0, inplace=True)

# Multidimensional Scaling
mds = MDS(n_components=2, dissimilarity='precomputed', random_state=1)
coords = mds.fit_transform(symmetrized_matrix)

# Hierarchical Clustering
linkage_matrix = linkage(coords, 'ward')

# Dendrogram
dendrogram(linkage_matrix, labels=symmetrized_matrix.index)
plt.title('D-Statistics Based Dendrogram')
plt.xlabel('Population')
plt.ylabel('Distance')
plt.show()


