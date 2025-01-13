# Interval Overlap Analysis Tool

## Overview
This tool is designed for genomic interval overlap analysis using interval trees and permutation testing. It calculates the number of overlapping bases between two sets of genomic intervals (`SetA` and `SetB`) and performs a permutation test to assess the statistical significance of the observed overlap.

---

## Features
- **Efficient Interval Tree Construction**: Organizes and queries intervals efficiently using a balanced interval tree.
- **Merge Overlapping Intervals**: Merges overlapping intervals in the input sets for accurate calculations.
- **Randomized Interval Permutations**: Randomizes intervals in `SetA` while maintaining interval lengths for permutation testing.
- **Permutation Test**: Calculates observed overlap and estimates p-values through repeated randomization.
- **File Parsing**: Supports parsing genomic interval files (BED format) and genome index files (`.fai` format).

---

## Usage

### Command-Line Interface
```bash
python interval_overlap_analysis.py SetA.bed SetB.bed genome.fa.fai [num_permutations]
