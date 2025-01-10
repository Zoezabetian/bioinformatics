# Protein Interaction Clustering and Visualization

## Overview
This script analyzes protein-protein interaction data to identify clusters of proteins based on their interaction scores. It uses non-negative matrix factorization (NMF) for dimensionality reduction, hierarchical clustering for grouping, and generates a heatmap with dendrograms for visualization.

---

## Features
1. **Load Data**:
   - Parses protein interaction data and optionally maps protein IDs to preferred names.
2. **Create Interaction Matrix**:
   - Constructs a symmetric matrix of interaction scores between proteins.
3. **Non-Negative Matrix Factorization (NMF)**:
   - Decomposes the interaction matrix to identify clusters based on protein similarity.
4. **Cluster Selection**:
   - Selects the smallest cluster of proteins for detailed analysis.
5. **Hierarchical Clustering**:
   - Uses linkage to group the most connected proteins in the selected cluster.
6. **Heatmap Visualization**:
   - Generates a heatmap with dendrograms to visualize patterns of protein interaction.

---

## Input Files
1. **Protein Links Data** (`protein_links.tsv`):
   - Contains interaction scores between pairs of proteins.
   - Format:
     ```
     protein1    protein2    combined_score
     ```
2. **Protein Info Data** (`protein_info.tsv`, optional):
   - Maps protein IDs to preferred names.
   - Format:
     ```
     #string_protein_id    preferred_name
     ```

---

## Output
1. **Heatmap**:
   - Saved as `Protein_Cluster_Heatmap.png`.
   - Displays interaction scores for top proteins in the selected cluster with hierarchical clustering dendrograms.

---

## Dependencies
- `numpy`: For numerical operations.
- `pandas`: For data manipulation and loading.
- `scikit-learn`: For NMF.
- `scipy`: For hierarchical clustering and distance computation.
- `seaborn`: For generating the heatmap.
- `matplotlib`: For visualization.

---

## Usage
Run the script in Python:
```bash
python protein_interaction_clustering.py <protein_links.tsv> [<protein_info.tsv>]
