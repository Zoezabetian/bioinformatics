import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage

# load the Dogs dataset
dogs_X = np.load('dogs_X.npy')
dogs_clades = np.load('dogs_clades.npy', allow_pickle=True)

# convert clades to numeric labels
unique_clades, clade_indices = np.unique(dogs_clades, return_inverse=True)

# function to compute clustering error
def compute_clustering_error(labels, true_labels):
    total_error = 0
    for label in np.unique(labels):
        cluster_labels = true_labels[labels == label]
        most_common = np.bincount(cluster_labels).argmax()
        misclassified = np.sum(cluster_labels != most_common)
        total_error += misclassified
    return total_error

# Perform hierarchical clustering
def hierarchical_clustering(X, true_labels, k):
    clustering = AgglomerativeClustering(n_clusters=k, linkage='average')
    labels = clustering.fit_predict(X)
    error = compute_clustering_error(labels, true_labels)
    return labels, error

# perform clustering and compute error
k = 30
labels, error = hierarchical_clustering(dogs_X, clade_indices, k)
print(f'Dogs Hierarchical Clustering Error (K={k}): {error}')

# generate the linkage matrix and create a dendrogram
linkage_matrix = linkage(dogs_X, method='average')
plt.figure(figsize=(12, 8))
dendro = dendrogram(
    linkage_matrix,
    truncate_mode='lastp',
    p=k,
    leaf_rotation=90,
    leaf_font_size=12,
    show_leaf_counts=True
)

# annotate with the most common clade
for i, (x, y) in enumerate(zip(dendro['icoord'], dendro['dcoord'])):
    x_pos = 0.5 * sum(x[1:3])
    y_pos = -0.5  # this is the vertical position where the label will be placed
    label = unique_clades[np.bincount(clade_indices[labels == i]).argmax()]
    plt.text(x_pos, y_pos, label, rotation=90, verticalalignment='top', horizontalalignment='center') # means 

plt.title('Dendrogram for Dogs Hierarchical Clustering (K=30)')
plt.xticks([])

plt.ylabel('Distance') 
plt.savefig('Dogs_dendrogram.png')
