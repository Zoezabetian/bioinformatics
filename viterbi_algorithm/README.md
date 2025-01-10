# Viterbi Algorithm for Identifying Inbred Regions in Genomic Data

## Overview
This script processes genomic data from a VCF file to identify regions of inbreeding using the Viterbi algorithm. It calculates the most likely sequence of inbred and outbred states for each individual across chromosomes and outputs the start and stop positions of inbred regions.

---

## Features
1. **VCF Parsing**:
   - Extracts genotype information, positions, and allele frequencies from a VCF file.
2. **Emission Probability Calculation**:
   - Computes the likelihood of observed genotypes for inbred and outbred states.
3. **Viterbi Algorithm**:
   - Identifies the most probable sequence of states (inbred/outbred) for each individual.
4. **Inbred Region Detection**:
   - Identifies and merges consecutive or adjacent inbred regions.

---

## Input
- **VCF File**: A variant call format (VCF) file containing genomic data.

---

## Output
- **Inbred Regions**:
  - Printed to the console in the format:
    ```
    individual    start    stop
    ```
  - Represents the genomic positions where inbreeding is most likely for each individual.

---

## Key Functions

### `parse_vcf(file_path)`
- Parses a VCF file and extracts genotype data, positions, and allele frequencies for each individual.

### `emission_prob(genotype, p, q, inbred_state)`
- Calculates the probability of observing a genotype in either the inbred or outbred state.

### `viterbi(genotypes, allele_freqs)`
- Implements the Viterbi algorithm to determine the most likely sequence of states (inbred/outbred) across positions.

### `main(vcf_file)`
- Processes the VCF file and outputs inbred regions for each individual.

---

## Dependencies
- `numpy`: For numerical computations.

---

## Usage
Run the script in Python:
```bash
python viterbi_algorithm.py synthetic_population.vcf
