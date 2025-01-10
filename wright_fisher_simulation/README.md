# Wright-Fisher Dynamics Simulation with Selection

## Overview
This script simulates the Wright-Fisher model with selection in a haploid population to study allele fixation or loss. The simulation is performed across multiple replicates, and the mean and variance of fixation and loss times are calculated.

---

## Features
1. **Wright-Fisher Simulation**:
   - Models allele frequency changes under selection until fixation (allele frequency = 1) or loss (allele frequency = 0).
2. **Fitness Adjustment**:
   - Incorporates relative fitness to adjust allele frequency probabilities.
3. **Simulation Outputs**:
   - Calculates the mean and variance of generations required for allele fixation or loss across replicates.

---

## Input Parameters
The script accepts the following command-line arguments:

- `--allele_freq`: Initial frequency of the allele (value between 0 and 1).
- `--pop_size`: Population size (integer).
- `--fitness`: Relative fitness of the allele.
- `--replicates`: Number of simulation replicates.

---

## Output
1. **Fixation Results**:
   - Mean and variance of generations required for the allele to fix.
2. **Loss Results**:
   - Mean and variance of generations required for the allele to be lost.

---

## Usage
Run the script in Python with the required arguments:
```bash
python wright_fisher_simulation.py --allele_freq <float> --pop_size <int> --fitness <float> --replicates <int>
