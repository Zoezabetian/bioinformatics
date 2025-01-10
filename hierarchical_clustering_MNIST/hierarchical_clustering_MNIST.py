# Zoe Zabetian Assignment 4 Tier 2
import numpy as np
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from scipy.cluster.hierarchy import dendrogram, linkage

# load the MNIST dataset subset
X = np.load('MNIST_X_subset.npy')
y = np.load('MNIST_y_subset.npy')

# function to compute clustering error
def compute_clustering_error(labels, y):
    total_error = 0
    unique_labels = np.unique(labels)  # find all unique clusters

    for cluster in unique_labels:
        # select the true labels of all points assigned to the current cluster
        cluster_labels = y[labels == cluster]
        if len(cluster_labels) > 0:
            # find the most common true label within the cluster
            most_common_label = np.bincount(cluster_labels).argmax()
            # count the number of misclassified points (i.e., points that don't match the most common label)
            misclassified_count = np.sum(cluster_labels != most_common_label)
            total_error += misclassified_count

    return total_error

# perform hierarchical clustering using agglomerative clustering
def hierarchical_clustering(X, y, k):
    # perform hierarchical clustering with average linkage
    clustering = AgglomerativeClustering(n_clusters=k, linkage='average')
    labels = clustering.fit_predict(X)  # fit and predict cluster labels

    # compute clustering error based on true labels
    error = compute_clustering_error(labels, y)

    return labels, error

# perform hierarchical clustering for k=10
k = 10
labels, error = hierarchical_clustering(X, y, k)

# save and display the clustering error
print(f'Hierarchical Clustering Error (K={k}): {error}')

# generate the linkage matrix for hierarchical clustering visualization
linkage_matrix = linkage(X, method='average')

# create a dendrogram with 10 terminal clusters (reduce number of leaves)
plt.figure(figsize=(12, 8))  # set figure size for clear visualization

dendro = dendrogram(
    linkage_matrix,
    truncate_mode='lastp',  # only show the last p merged clusters
    p=k,  # number of clusters to show
    leaf_rotation=90,
    leaf_font_size=10,
    show_leaf_counts=True
)

# Annotate with the most common label for each cluster
for i, (x, y_coords) in enumerate(zip(dendro['icoord'], dendro['dcoord'])):
    x_pos = 0.5 * (x[1] + x[2])  # middle of the cluster
    y_pos = y_coords[1]  # vertical position of the cluster
    label = np.bincount(y[labels == i]).argmax()
    plt.text(
        x=x_pos,
        y=y_pos,
        s=str(label),
        horizontalalignment='center',
        verticalalignment='top',
        fontsize=10,
        color='red'
    )

# set title and axis labels
plt.title('Dendrogram for MNIST Hierarchical Clustering (K=10)')

plt.ylabel('Distance')

# save the dendrogram as a png file
plt.xticks([])
plt.savefig('MNIST_dendrogram.png')

# save explanation to a text file
with open('MNIST_paragraph.txt', 'w') as f:
    f.write("The dendrogram shows hierarchical clustering of the MNIST digits, with the clusters representing the closest groups based on pixel similarity. I couldn't figure out how to get the labels on the terminal nodes, but from what I do notice from the labels at the roots, the numbers that are more similar in shape are closer together, such as the 1 and 7, which are next to each other, as well as the 4 and 9 (for example). This shows that the clustering algorithm captures visual similarity between the digits.")