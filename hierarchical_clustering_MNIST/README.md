# Hierarchical Clustering on MNIST Subset

## Overview
This script performs hierarchical clustering on a subset of the MNIST dataset using agglomerative clustering. It visualizes the clustering process with a dendrogram and calculates clustering error based on the true digit labels.

---

## Features
- **Hierarchical Clustering**:
  - Implements agglomerative clustering with average linkage for `k=10` clusters.
- **Dendrogram Visualization**:
  - Generates a dendrogram to visualize hierarchical relationships between clusters.
- **Clustering Error Calculation**:
  - Calculates the clustering error based on the most common label within each cluster.
- **Annotated Dendrogram**:
  - Displays the most common label for each cluster at key points in the dendrogram.
- **Explanatory Text**:
  - Saves an explanation of clustering results to a text file.

---

## Input Data
- **`MNIST_X_subset.npy`**: A NumPy file containing flattened 28x28 pixel data for 6000 images (each image is a 784-element array).
- **`MNIST_y_subset.npy`**: A NumPy file containing true digit labels (0-9) for the corresponding images.

---

## Output
1. **Dendrogram**:
   - Saved as `MNIST_dendrogram.png`, showing hierarchical clustering relationships.
2. **Clustering Error**:
   - Printed in the console:
     ```
     Hierarchical Clustering Error (K=10): ...
     ```
3. **Explanatory Text**:
   - Saved in `MNIST_paragraph.txt` with insights on clustering results and visual similarities between digits.

---

## Key Functions

### `hierarchical_clustering(X, y, k)`
- **Input**:
  - `X`: The input dataset (images as flattened arrays).
  - `y`: True digit labels.
  - `k`: Number of clusters.
- **Output**:
  - `labels`: Predicted cluster labels for each data point.
  - `error`: Total clustering error based on true labels.
- **Description**:
  - Performs agglomerative clustering with average linkage.
  - Assigns each data point to a cluster and calculates misclassified points.

### `compute_clustering_error(labels, y)`
- **Input**:
  - `labels`: Predicted cluster labels.
  - `y`: True digit labels.
- **Output**:
  - `total_error`: Number of misclassified points across all clusters.
- **Description**:
  - Identifies the most common true label in each cluster and counts misclassified points.

---

## Workflow
1. **Data Loading**:
   - Load `MNIST_X_subset.npy` and `MNIST_y_subset.npy`.
2. **Hierarchical Clustering**:
   - Perform clustering for `k=10`.
   - Compute clustering error.
3. **Dendrogram Visualization**:
   - Generate and save a dendrogram with hierarchical relationships between clusters.
4. **Explanatory Text**:
   - Save observations about clustering results to `MNIST_paragraph.txt`.

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
python hierarchical_clustering_MNIST.py
