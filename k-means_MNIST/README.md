# K-Means Clustering on MNIST Subset

## Overview
This script performs k-means clustering on a subset of the MNIST dataset, which contains images of handwritten digits. The goal is to cluster the images and evaluate the clustering error based on true digit labels.

---

## Features
- Loads and visualizes MNIST image data.
- Implements k-means clustering with `k=10` and `k=11`.
- Calculates clustering error based on misclassified points in each cluster.
- Visualizes cluster centroids as grayscale images.
- Saves centroid images for each value of `k` as PNG files.

---

## Input Data
- **`MNIST_X_subset.npy`**: A NumPy file containing flattened 28x28 pixel data for 6000 images (each image is a 784-element array).
- **`MNIST_y_subset.npy`**: A NumPy file containing true digit labels (0-9) for the corresponding images.

---

## Output
1. **Cluster Centroid Images**:
   - For each `k`, the script saves centroid images as:
     - `centroids_k10.png`
     - `centroids_k11.png`
2. **Clustering Error**:
   - Prints clustering error for each value of `k`:
     ```
     K=10 Error=...
     K=11 Error=...
     ```

---

## Key Functions

### `kmeans_clustering(X, y, k)`
- **Input**:
  - `X`: The input dataset (images as flattened arrays).
  - `y`: True digit labels.
  - `k`: Number of clusters.
- **Output**:
  - `centroids`: The mean positions of clusters.
  - `total_error`: The number of misclassified points across all clusters.
- **Description**:
  - Fits a k-means model to the dataset.
  - Calculates the most common label in each cluster.
  - Counts and sums misclassified points for all clusters.

---

## Workflow
1. **Data Loading**:
   - Load `MNIST_X_subset.npy` and `MNIST_y_subset.npy`.
2. **Data Visualization**:
   - Display the first example from the dataset.
3. **K-Means Clustering**:
   - Perform clustering for `k=10` and `k=11`.
   - Reshape and save centroids as images.
4. **Error Calculation**:
   - Compute and print clustering error for each `k`.

---

## Dependencies
- `numpy`: For numerical operations and data handling.
- `matplotlib`: For visualizing images and centroids.
- `scikit-learn`: For implementing k-means clustering.

---

## Usage
Run the script in Python:
```bash
python k-means_MNIST.py
