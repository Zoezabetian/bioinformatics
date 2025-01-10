# Fibonacci Sequence Calculator

## Overview
This script calculates the `n`-th term of a generalized Fibonacci-like sequence using a closed-form formula, avoiding recursion or iteration. It leverages the mathematical properties of eigenvalues and the golden ratio (φ) for efficient computation.

---

## Features
- **Closed-Form Calculation**:
  - Computes the `n`-th term without generating all preceding terms.
  - Utilizes the golden ratio (φ) and its conjugate (ψ) to derive a direct solution.
- **Supports Custom Starting Values**:
  - The sequence is initialized with user-defined terms `a` and `b`.

---

## Input
- **Command-Line Arguments**:
  - `a`: First term of the sequence.
  - `b`: Second term of the sequence.
  - `n`: Term index to calculate.

---

## Output
- Prints the `n`-th term of the sequence.

---

## Formula
The script calculates the `n`-th term as: Fn = c1 * (φ^(n-1)) + c2 * (ψ^(n-1))
Where:
- `φ = (1 + √5) / 2`: The golden ratio.
- `ψ = (1 - √5) / 2`: The conjugate of φ.
- `c1` and `c2`: Coefficients derived from the initial values `a` and `b`.

---

## Usage
Run the script in Python with the following command:
```bash
python fibonacci_calculator.py <a> <b> <n>

