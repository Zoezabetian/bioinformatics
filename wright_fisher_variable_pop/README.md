# Variable Population Wright-Fisher Simulation with Selection

## Overview
This script simulates the Wright-Fisher model with selection in a haploid population under variable population sizes over generations. The population size changes are provided through a TSV file, enabling flexible simulation dynamics. The script computes the number of generations required for an allele to fix (frequency = 1) or be lost (frequency = 0) across multiple replicates.

---

## Features
1. **Variable Population Sizes**:
   - Allows population size to change across generations based on an input schedule.
2. **Fitness Adjustment**:
   - Incorporates relative fitness to adjust allele inheritance probabilities.
3. **Simulation Outputs**:
   - Calculates the mean and variance of fixation and loss times across replicates.

---

## Input Parameters
The script accepts the following command-line arguments:

- `--allele_freq`: Initial frequency of the allele (value between 0 and 1).
- `--fitness`: Relative fitness of the allele.
- `--replicates`: Number of simulation replicates.
- `--pop_size_file`: Path to the TSV file specifying the population size schedule.

### Population Size File Format
The TSV file should have two columns: `generation` and `popsize` (tab-separated). The `generation` column specifies the generation where the population size changes, and the `popsize` column specifies the corresponding population size.

## Usage
Run the script in Python with the required arguments:
```bash
python wright_fisher_variable_pop.py --allele_freq <float> --fitness <float> --replicates <int> --pop_size_file <path>
