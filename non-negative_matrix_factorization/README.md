# Non-Negative Matrix Factorization (NMF) on Dog SNP Data

## Overview
This script analyzes SNP data from dogs to identify genetic patterns using Non-Negative Matrix Factorization (NMF). It decomposes the dataset into components that represent shared genetic traits among dogs and visualizes the proportions of these components for each individual.

---

## Features
1. **Load Data**:
   - Loads SNP data and corresponding clade labels for dogs.
2. **Apply NMF**:
   - Decomposes the SNP dataset into two matrices:
     - **Matrix W**: Represents each dog as a combination of components.
     - **Matrix H**: Represents how SNP features contribute to each component.
3. **Normalize Components**:
   - Converts each dog's representation into proportions across components for easier interpretation.
4. **Sort by Dominant Cluster**:
   - Groups and sorts dogs by their dominant component and its proportion.
5. **Stacked Plot Visualization**:
   - Generates a stacked bar plot showing the proportions of components for each dog.
6. **Dominant Component Identification**:
   - Identifies the dominant components in Basenji and Wolf samples.

---

## Input Files
1. **Dog SNP Data** (`dog_X.npy`):
   - Each row represents a dog, and each column represents an SNP feature.
2. **Clade Labels** (`dog_clades.npy`):
   - Labels identifying the clade of each dog.

---

## Output
1. **Visualization**:
   - Stacked bar plot saved as `NMF_Dogs.png`:
     - X-axis: Individual dogs, sorted by dominant cluster.
     - Y-axis: Proportion of components.
2. **Console Output** (Optional):
   - Dominant components for Basenji and Wolf samples.

---

## Dependencies
- `numpy`: For numerical operations.
- `pandas`: For data manipulation.
- `scikit-learn`: For NMF decomposition.
- `matplotlib`: For visualization.

---

## Usage
Run the script in Python:
```bash
python non-negative_matrix_factorization.py <dog_X.npy> <dog_clades.npy>
