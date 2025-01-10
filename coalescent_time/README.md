# Coalescent Time Simulation

## Overview
This script simulates the time to the eighth coalescent event in a Wright-Fisher population model. It computes the average time and variance of the specified coalescent event across multiple replicates.

## Features
- Simulates coalescent events with exponential waiting times.
- Models lineage reduction over time due to coalescent events.
- Configurable population size, sample size, and number of replicates.

## Input Parameters
- `--pop_size`: Effective population size (integer).
- `--sample_size`: Initial number of lineages (integer).
- `--replicates`: Number of simulation replicates (integer).

## Output
The mean and variance of the time required to reach the eighth coalescent event.

## Usage
Run the script with the required arguments:
```bash
python coalescent_time.py --pop_size <int> --sample_size <int> --replicates <int>
