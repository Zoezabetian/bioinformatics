# 3D Molecular Visualization

## Overview
This script visualizes a molecular structure in 3D using MDS coordinates and an atomic distance matrix. Atoms are color-coded and sized based on CPK conventions, and bonds are drawn between atoms based on their distances.

---

## Features
- **3D Molecular Visualization**:
  - Uses 3D scatter plotting to represent atoms.
  - Color-codes atoms based on their element using the CPK color convention.
  - Sizes atoms proportionally to their atomic radii.
- **Bond Representation**:
  - Draws bonds between atoms with distances less than 1.6 units.
- **Legend**:
  - Includes a legend indicating atom types and their corresponding colors.

---

## Input Data
- **`molecule_coordinates.csv`**:
  - Contains MDS coordinates (`x`, `y`, `z`) and atom types (`element`).
  - Columns:
    ```
    x, y, z, element
    ```
- **`molecule_distances.tsv`**:
  - Contains a molecular distance matrix, where each value represents the distance between a pair of atoms.

---

## Output
- **3D Visualization**:
  - Saved as `molecule_3D_plot.png`, a 3D scatter plot with atoms and bonds.

---

## Key Details

### CPK Colors and Atomic Sizes
- **Colors**:
  - Hydrogen: White
  - Carbon: Black
  - Oxygen: Red
- **Sizes**:
  - Hydrogen: 50
  - Carbon: 100
  - Oxygen: 140

### Bonds
- Bonds are drawn between atoms with distances less than 1.6 units.

---

## Dependencies
- `numpy`: For numerical operations and data handling.
- `pandas`: For loading and processing input data.
- `matplotlib`: For 3D visualization.

---

## Usage
Run the script in Python:
```bash
python 3D_molecular_visualization.py
