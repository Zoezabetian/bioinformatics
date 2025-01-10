# Hierarchical Clustering on Dogs Dataset

## Overview
This script performs hierarchical clustering on a dataset of dog breeds (`dogs_X.npy`), classifies them into clades (`dogs_clades.npy`), and evaluates clustering performance using a dendrogram visualization.

---

## Features
- **Hierarchical Clustering**:
  - Implements agglomerative clustering with average linkage for `k=30` clusters.
- **Clustering Error Calculation**:
  - Computes the error based on misclassified clade labels.
- **Dendrogram Visualization**:
  - Generates a dendrogram to illustrate hierarchical clustering relationships.
- **Clade Annotation**:
  - Annotates the dendrogram with the most common clade label for each cluster.

---

## Input Data
- **`dogs_X.npy`**: A NumPy file containing feature data for dog breeds.
- **`dogs_clades.npy`**: A NumPy file containing clade labels for corresponding dog breeds.

---

## Output
1. **Dendrogram**:
   - Saved as `Dogs_dendrogram.png`, showing hierarchical relationships between clusters and their dominant clades.
2. **Clustering Error**:
   - Printed in the console:
     ```
     Dogs Hierarchical Clustering Error (K=30): ...
     ```

---

## Key Functions

### `compute_clustering_error(labels, true_labels)`
- **Input**:
  - `labels`: Predicted cluster labels.
  - `true_labels`: True clade labels.
- **Output**:
  - `total_error`: Number of misclassified points across all clusters.
- **Description**:
  - For each cluster, identifies the most common clade and counts misclassified points.

### `hierarchical_clustering(X, true_labels, k)`
- **Input**:
  - `X`: Feature data for the dataset.
  - `true_labels`: True clade labels.
  - `k`: Number of clusters.
- **Output**:
  - `labels`: Predicted cluster labels.
  - `error`: Total clustering error.
- **Description**:
  - Performs agglomerative clustering and computes the clustering error.

---

## Workflow
1. **Data Loading**:
   - Load `dogs_X.npy` and `dogs_clades.npy`.
   - Convert clade labels to numeric values for clustering.
2. **Hierarchical Clustering**:
   - Perform clustering for `k=30`.
   - Compute clustering error based on true clade labels.
3. **Dendrogram Visualization**:
   - Generate and save a dendrogram with hierarchical relationships between clusters.
   - Annotate dendrogram clusters with their dominant clades.

---

## Dependencies
- `numpy`: For numerical operations and data handling.
- `matplotlib`: For generating dendrogram visualizations.
- `scikit-learn`: For implementing agglomerative clustering.
- `scipy`: For generating linkage matrices.

---

## Usage
Run the script in Python:
```bash
python hierarchical_clustering_dogs.py
