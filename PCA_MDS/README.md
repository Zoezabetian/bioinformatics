# PCA and MDS Analysis on Multiple Datasets

## Overview
This script performs dimensionality reduction using PCA and MDS on three datasets: MNIST, Dogs SNP data, and a molecule distance matrix. It visualizes the results and saves outputs including transformed data, reconstructed examples, and scatter plots.

---

## Features

### **Part 1: PCA on MNIST Dataset**
- Reduces MNIST digit data from 784 dimensions to 2 dimensions.
- Visualizes the 2D PCA results with digit labels.
- Reconstructs a sample image using the first 2 principal components.
- Saves:
  - Original and reconstructed images (`MNIST_original.png`, `MNIST_reconstructed_2PC.png`).
  - Reconstructed image from a specific 2D coordinate (`MNIST_reconstructed_1_from_coord.png`).
  - 2D PCA scatter plot (`MNIST_PCA_2D.png`).

### **Part 2: PCA on Dogs SNP Dataset**
- Reduces SNP data for dog breeds from 784 dimensions to 2 dimensions.
- Visualizes clades in a 2D PCA scatter plot, color-coded by clade labels.
- Saves:
  - 2D PCA scatter plot (`Dogs_PCA_2D.png`).

### **Part 3: MDS on Molecule Distance Matrix**
- Applies MDS to a molecular distance matrix to obtain 3D coordinates.
- Saves:
  - 3D molecule coordinates to a CSV (`molecule_coordinates.csv`).

---

## Input Data

### MNIST Dataset
- **`MNIST_X_subset.npy`**: A NumPy file containing 784-dimensional feature vectors for MNIST digits.
- **`MNIST_y_subset.npy`**: A NumPy file containing digit labels.

### Dogs SNP Dataset
- **`dogs_X.npy`**: A NumPy file containing SNP feature vectors for dog breeds.
- **`dogs_clades.npy`**: A NumPy file containing clade labels for dog breeds.

### Molecule Distance Matrix
- **`molecule_distances.tsv`**: A tab-separated file containing a molecular distance matrix (NMR data).

---

## Outputs

### Part 1: MNIST PCA
- `MNIST_PCA_2D.png`: 2D scatter plot of PCA-transformed MNIST data.
- `MNIST_original.png`: Original MNIST image (28x28 pixels).
- `MNIST_reconstructed_2PC.png`: Reconstructed image using the first 2 principal components.
- `MNIST_reconstructed_1_from_coord.png`: Reconstructed image from specific 2D coordinates.

### Part 2: Dogs SNP PCA
- `Dogs_PCA_2D.png`: 2D scatter plot of PCA-transformed SNP data, color-coded by clade.

### Part 3: Molecule MDS
- `molecule_coordinates.csv`: CSV file containing 3D coordinates of molecules from MDS.

---

## Dependencies
- `numpy`: For numerical operations.
- `pandas`: For data handling and I/O.
- `matplotlib`: For data visualization.
- `scikit-learn`: For PCA and MDS.

---

## Usage
Run the script in Python:
```bash
python PCA_MDS.py
