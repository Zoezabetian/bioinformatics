# Gene Expression Analysis Tools

## Overview
This folder contains two Python scripts for analyzing differential gene expression under control and treatment conditions.

---

## Scripts

`differential_gene_expression.py`
- Purpose: Identifies differentially expressed genes based on normalized counts and log₂ fold changes.
- Output: `output1.txt` with columns:
  gene  mean_normalized_control_expression  median_normalized_control_expression
  mean_normalized_treatment_expression  median_normalized_treatment_expression  logFoldChange
- Usage:
  python differential_gene_expression.py path/to/control_dir path/to/treatment_dir

`gene_expression_statistical_testing.py`
- Purpose: Extends analysis with Mann-Whitney U test for statistical significance of median expression differences.
- Output: `output2.txt` with columns:
  gene  mean_normalized_control_expression  median_normalized_control_expression
  mean_normalized_treatment_expression  median_normalized_treatment_expression  logFoldChange  p_value
- Usage:
  python gene_expression_statistical_testing.py path/to/control_dir path/to/treatment_dir

---

## Input Data
- Format: CSV files with:
  gene_name,expression_value
- Directories:
  control_dir: CSVs for the control group.
  treatment_dir: CSVs for the treatment group.
