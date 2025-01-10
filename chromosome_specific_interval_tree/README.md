# Genomic Interval Overlap Analysis Tool

## Overview
This tool is designed to calculate the overlap between two sets of genomic intervals (`SetA` and `SetB`) and perform a permutation test to evaluate the statistical significance of the observed overlap. It uses interval trees for efficient overlap calculation and supports randomized interval generation for permutation testing.

---

## Features
- **Efficient Interval Tree Construction**: Organizes intervals into a balanced interval tree for fast overlap queries.
- **Interval Merging**: Merges overlapping intervals within a set for accurate calculations.
- **Chromosome-Specific Processing**: Handles genomic intervals on a per-chromosome basis.
- **Permutation Testing**: Generates randomized intervals to compute p-values for observed overlaps.
- **Parallel Chromosome Analysis**: Separates interval processing by chromosome for efficient computation.

---

## Usage

### Command-Line Interface
```bash
python chromosome_specific_interval_tree.py path/to/SetA.bed path/to/SetB.bed path/to/genome.fa.fai [num_permutations]
