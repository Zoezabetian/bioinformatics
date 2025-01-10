import numpy as np
import pandas as pd
from sklearn.decomposition import NMF
from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import warnings

warnings.simplefilter("ignore", category=UserWarning)  # ignore seaborn clustermap warning

def load_data(filepath_links, filepath_info=None):
    """Load the protein-protein interaction dataset and optionally protein info."""
    data = pd.read_csv(filepath_links, sep='\t')  # load protein-protein interaction dataset
    
    # load protein info only if filepath_info is provided
    if filepath_info:
        info = pd.read_csv(filepath_info, sep='\t')  # load protein info dataset
        id_to_name = dict(zip(info['#string_protein_id'], info['preferred_name']))  # map protein IDs to names
    else:
        id_to_name = None  # no ID mapping if info file is not provided
    return data, id_to_name


def create_interaction_matrix(data):
    """constructs a symmetric matrix where each protein is linked by interaction scores, capturing the connections between proteins."""
    proteins = sorted(set(data['protein1']).union(set(data['protein2'])))  # get unique proteins and sort them
    protein_idx = {protein: i for i, protein in enumerate(proteins)}  # map protein names to indices
    n_proteins = len(proteins)  # get the number of proteins

    # initialize a matrix of zeros, representing no interaction until populated
    interaction_matrix = np.zeros((n_proteins, n_proteins)) 

    for _, row in data.iterrows():
        i, j = protein_idx[row['protein1']], protein_idx[row['protein2']]  # get the indices of the two proteins
        interaction_matrix[i, j] = row['combined_score']  # set the interaction scores
        interaction_matrix[j, i] = row['combined_score']  
        # value at i,j and j,i should be the same as the matrix is symmetric
        # we want symmetry because if one protein interacts with another, the other interacts with the first
    return interaction_matrix, proteins


def apply_nmf(interaction_matrix):
    """Apply NMF to decompose the interaction matrix. Get matrix W, with dimensions (n_proteins, n_components), where 
    each row (protein) contains values across 10 components, indicating its membership strength in each component
    """
    #  dimensionality reduction technique that factorizes a matrix into two lower-rank matrices: W (basis matrix) and H (coefficient matrix)
    nmf = NMF(n_components=10, init='random', random_state=42)  # initialize NMF model, 10 clusters
    W = nmf.fit_transform(interaction_matrix)  # fit the model and transform the matrix
    return W # Only need W because you get the strength of each protein’s association with each component (where each component can be 
             # thought of as a type of cluster), and cluster assignments for each protein, based on the highest value in each row of W
            

def assign_clusters(W):
    """Each protein is assigned to a cluster based on the highest value in W, grouping proteins by similarity
    """
    cluster_assignments = np.argmax(W, axis=1)  # assigns each protein to the cluster (component) where it has the highest membership value
    cluster_counts = np.bincount(cluster_assignments)  # count the number of proteins in each cluster
    smallest_cluster = np.argmin(cluster_counts)  # find the smallest cluster (fewest proteins)
    selected_cluster_proteins = np.where(cluster_assignments == smallest_cluster)[0]  # get the indices of proteins in the smallest cluster
    return selected_cluster_proteins


def calculate_linkage(cluster_matrix, top_n=50):
    """using the selected cluster, proteins with the highest total interactions are chosen, and a linkage matrix is calculated to perform hierarchical clustering
    Linkage = method used in hierarchical clustering to determine how clusters are formed by combining individual data points or existing clusters based on a chosen distance measure. 
    In hierarchical clustering, linkage defines the rules for calculating the "distance" between clusters at each step, helping to build a tree-like structure called a dendrogram.
    """
    interaction_sums = cluster_matrix.sum(axis=1)  # sum the interactions for each protein
    top_indices = interaction_sums.argsort()[-top_n:]  # proteins with the highest interaction sums are chosen as the most "connected" proteins. -top_n: get the top n indices
    top_cluster_matrix = cluster_matrix[np.ix_(top_indices, top_indices)]  # subset the matrix to the top proteins. np.ix_ is used to index the matrix
    
    np.fill_diagonal(top_cluster_matrix, 0)  # set the diagonal to zero (no self-interactions)

    condensed_matrix = squareform(top_cluster_matrix)  # squareform converts the symmetric interaction matrix into a condensed distance format required by linkage
    linkage_matrix = linkage(condensed_matrix, method='average')  # performs average linkage clustering, where each cluster’s distance is the average distance between members
    return top_cluster_matrix, top_indices, linkage_matrix

def plot_heatmap(top_cluster_matrix, top_protein_names, linkage_matrix):
    """Plot a heatmap with dendrogram for the selected cluster of proteins.
    """
    g = sns.clustermap(
        top_cluster_matrix,
        cmap='viridis',
        row_cluster=True,        # enables hierarchical clustering on rows, resulting in a row dendrogram
        col_cluster=True,        # enables hierarchical clustering on columns, resulting in a column dendrogram
        row_linkage=linkage_matrix,  # specifies the linkage matrix for clustering the rows
        col_linkage=linkage_matrix,  # specifies the linkage matrix for clustering the columns
        xticklabels=top_protein_names,
        yticklabels=top_protein_names,
        figsize=(18, 18)
    )
    g.ax_heatmap.set_xticklabels(g.ax_heatmap.get_xmajorticklabels(), fontsize=8, rotation=45)
    g.ax_heatmap.set_yticklabels(g.ax_heatmap.get_ymajorticklabels(), fontsize=8)
    plt.savefig("Protein_Cluster_Heatmap.png")


def main(filepath_links, filepath_info=None):
    data, id_to_name = load_data(filepath_links, filepath_info)  # load the data and map protein IDs to names if available
    interaction_matrix, proteins = create_interaction_matrix(data)  # constructs a symmetric matrix where each protein is linked by interaction scores, capturing the connections between proteins.
    W = apply_nmf(interaction_matrix)  # reduce data into two matrices—one showing each protein’s combination of components (W) and the other defining the components themselves (H).
    selected_cluster = assign_clusters(W)  # Each protein is assigned to a cluster based on the highest value in W, grouping proteins by similarity
    selected_cluster_matrix = interaction_matrix[np.ix_(selected_cluster, selected_cluster)]  # the smallest cluster is selected for further analysis
    top_cluster_matrix, top_indices, linkage_matrix = calculate_linkage(selected_cluster_matrix)  #  using the selected cluster, proteins with the highest total interactions are chosen, and a linkage matrix is calculated to perform hierarchical clustering.
    
    # use ID mapping if available, otherwise keep original IDs
    top_protein_names = [
        id_to_name.get(proteins[selected_cluster[i]], proteins[selected_cluster[i]]) if id_to_name 
        else proteins[selected_cluster[i]]
        for i in top_indices
    ]
    
    plot_heatmap(top_cluster_matrix, top_protein_names, linkage_matrix)  # heatmap with clustering dendrograms is generated, visualizing patterns of interaction among proteins.


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python submission.py <protein_links.tsv> [<protein_info.tsv>]")
    elif len(sys.argv) == 2:
        main(sys.argv[1])  # only protein links provided
    else:
        main(sys.argv[1], sys.argv[2])  # both protein links and info files provided


# Why don't we use Matrix H? It could help if we wanted to interpret the nature of each component (e.g., understanding which features define each cluster). 
# However, this is not required for simply clustering and visualizing protein relationships.