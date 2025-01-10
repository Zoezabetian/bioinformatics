# Transition Rate Optimization Using Forward Algorithm and Amoeba Optimization

## Overview
This script calculates transition rates between inbred and outbred states using genotype data from a VCF file. The Forward Algorithm is used to compute the likelihood of the observed data, and the Amoeba optimization algorithm estimates the optimal transition rates.

---

## Features
- **VCF Parsing**:
  - Reads genotype data from a VCF file.
- **Emission Probabilities**:
  - Computes the likelihood of genotypes given inbred and outbred states.
- **Forward Algorithm**:
  - Computes the total likelihood of observed genotypes for given transition rates.
- **Amoeba Optimization**:
  - Optimizes transition rates to maximize the likelihood of observed data.

---

## Input
- **VCF File**: A variant call format file containing genotype data.
- **Initial Parameters**:
  - `P(transition outbred → inbred)`: Default initial value is `1 / (4 * 10^6)`.
  - `P(transition inbred → outbred)`: Default initial value is `1 / (1.5 * 10^6)`.

---

## Output
- Optimized transition rates:
  - `P(transition outbred → inbred)`
  - `P(transition inbred → outbred)`

---

## Usage
Run the script in Python:
```bash
python amoeba_optimization.py synthetic_population.vcf
