# Zoe's Bioinformatics Projects

## Overview
This repository contains a collection of computational biology and machine learning projects developed for academic, research, and fun purposes. The projects span a wide range of topics, including clustering, machine learning, population genetics, and molecular visualization. Each project includes scripts, data analysis workflows, and methods to address specific biological problems or computational challenges.

---

## Projects:

### 1. **3D Molecular Visualization**
- **Description:** Visualize molecular structures in 3D using a provided distance matrix.
- **Key Features:**
  - Uses CPK coloring for atoms.
  - Generates bonds based on distances between atoms.
  - Saves a 3D molecular visualization as an image.

---

### 2. **PCA and MDS**
- **Description:** Dimensionality reduction on MNIST and molecular data using PCA and MDS.
- **Key Features:**
  - Visualizes MNIST digits in 2D.
  - Reconstructs MNIST images using reduced dimensions.
  - Applies MDS on molecular distances.

---

### 3. **XGBoost Breast Cancer Prediction**
- **Description:** Predict breast cancer using recursive feature elimination and hyperparameter tuning with XGBoost.
- **Key Features:**
  - Recursive feature elimination for feature selection.
  - Cross-validation with hyperparameter grid search.
  - Predicts probabilities for test samples.

---

### 4. **Amoeba Optimization**
- **Description:** Uses the Amoeba optimization algorithm to adjust parameters for transition probabilities in a bioinformatics model.
- **Key Features:**
  - Applies forward algorithms for likelihood computation.
  - Optimizes probabilities for inbred and outbred states.

---

### 5. **Bioinformatics Practice Notebook**
- **Description:** A Jupyter Notebook for exploring various bioinformatics problems from Rosalind.

---

### 6. **Chromosome-Specific Interval Tree**
- **Description:** Builds an interval tree for chromosome-specific analyses and finds overlapping intervals.
- **Key Features:**
  - Efficient overlap calculation.
  - Supports querying chromosome intervals.

---

### 7. **Coalescent Time Simulation**
- **Description:** Simulates the time to coalescent events in populations.
- **Key Features:**
  - Models population dynamics for a Wright-Fisher population.
  - Calculates the time to specific coalescent events.

---

### 8. **Differential Gene Expression Analysis**
- **Description:** Analyzes differential gene expression using log fold changes and statistical tests.
- **Key Features:**
  - Computes normalized expression and fold changes.
  - Sorts genes by significance.

---

### 9. **DNA as Data Storage: Encoding and Decoding**
- **Description:** Encodes binary files into tagged DNA sequences and decodes them back to their original format.
- **Key Features:**
  - Inserts location-specific tags to improve data segmentation.
  - Ensures sequence consistency with non-coding filler nucleotides.
  - Demonstrates robust encoding and decoding for data integrity.

---

### 10. **Fibonacci Calculator**
- **Description:** Computes Fibonacci numbers using a closed-form solution.
- **Key Features:**
  - Implements matrix-based representation for efficient calculation.

---

### 11. **Hierarchical Clustering (MNIST)**
- **Description:** Applies hierarchical clustering on MNIST digits and visualizes the dendrogram.
- **Key Features:**
  - Clusters digits based on visual similarity.
  - Saves a dendrogram as an image.
- **Output:** `MNIST_dendrogram.png`.

---

### 12. **Hierarchical Clustering (Dogs)**
- **Description:** Clusters dog breeds based on genetic similarities.
- **Key Features:**
  - Creates a dendrogram for hierarchical clustering.
  - Saves dendrogram for genetic insights.
- **Output:** `Dogs_dendrogram.png`.

---

### 13. **Interval Overlap Analysis**
- **Description:** Analyzes overlaps between genomic intervals using permutation testing.
- **Key Features:**
  - Randomizes interval positions.
  - Calculates observed and permuted overlaps.

---

### 14. **K-Means Clustering (MNIST)**
- **Description:** Applies K-Means clustering to MNIST digit data.
- **Key Features:**
  - Visualizes centroids as digit-like images.
  - Calculates clustering errors.

---

### 15. **Non-Negative Matrix Factorization (NMF)**
- **Description:** Decomposes SNP data to reveal shared patterns among dogs.
- **Key Features:**
  - Identifies dominant clusters for each dog.
  - Visualizes cluster proportions as a stacked plot.

---

### 16. **Protein Interaction Clustering**
- **Description:** Clusters proteins based on interaction scores using hierarchical clustering.
- **Key Features:**
  - Generates a heatmap with dendrograms for clustering visualization.
  - Identifies clusters among the top interacting proteins.

---

### 17. **Viterbi Algorithm**
- **Description:** Applies the Viterbi algorithm to predict hidden states in genetic sequences.
- **Key Features:**
  - Parses VCF files for genotypes and frequencies.
  - Predicts inbred and outbred regions across chromosomes.

---

### 18. **Wright-Fisher Simulation**
- **Description:** Simulates allele fixation and loss under Wright-Fisher dynamics with constant population sizes.
- **Key Features:**
  - Models the effect of selection on allele frequencies.

---

### 19. **Wright-Fisher Simulation with Variable Population Sizes**
- **Description:** Extends the Wright-Fisher model to include dynamic population sizes.
- **Key Features:**
  - Incorporates population schedules from TSV files.
  - Tracks fixation and loss events over generations.
