# Optimization Algorithm ——IGKSO ——improved gkso

## Description
This repository provides a comprehensive framework for benchmarking and analyzing the performance of intelligent optimization algorithms. It includes implementations of statistical tests (Friedman test, Wilcoxon rank-sum test), algorithm code (IGKSOne), standard test suites (CEC2017, CEC2019), and tools for visualizing iterative convergence behavior. The codebase is designed for researchers in evolutionary computation and swarm intelligence to conduct reproducible experiments and comparative studies.

## Dataset Information
- **CEC2017 Test Suite**: 29 benchmark functions (F1–F30 excluding F2) with configurable dimensions (30D). 
- **CEC2019 Test Suite**: 10 benchmark functions (F1–F10 ). 

## Code Information
The repository contains the following components:
- **Fridman test code** – Implements Friedman test for multiple algorithm comparison.
- **IGKSOne code** – Implementation of Improved Genghis Khan Shark Optimization Algorithm.
- **Wilcoxon rank-sum test results** – Implements Wilcoxon rank-sum test code for multiple algorithm comparison.
- **cec2017 test suite code** – Integration of CEC2017 benchmark functions.
- **cec2019 test suite code** – Integration of CEC2019 benchmark functions.
- **iterative graph code** – iterative graph code.

## Usage Instructions

1. **Run Algorithms**  
   Use `IGKSOne code` to execute the optimization algorithm on selected CEC functions from `cec2017 test suite code` or `cec2019 test suite code`.
   Select main.m in the file to execute code.
   Save performance metrics (best fitness per iteration, mean and std) for statistical analysis.
2. **Perform Statistical Tests**  
   - Execute `Fridman test code` for overall ranking of multiple algorithms.
      Select main.m in the file to execute code. 
   - Review `Wilcoxon rank-sum test results` for pairwise significance testing.
       Select zhihejianyan.m in the file to execute code.
3. **Visualize Results**  
   Use `iterative graph code` to plot convergence curves and diversity trends.
   Select f1.m in the file to execute code.


## Requirements
- **Programming Language**: MATLAB

## Methodology
1. **Initialization**: Configure algorithm parameters and test function settings (dimension, run count).
2. **Execution**: Run optimization algorithm on each test function for a predefined number of iterations/runs.
3. **Data Recording**: Capture best-so-far fitness,mean,std, convergence curves, and other statistics.
4. **Statistical Analysis**:
   - Apply Friedman test for multiple algorithm ranking.
   - Apply Wilcoxon rank-sum test for pairwise comparisons.
5. **Visualization**: Generate plots to illustrate convergence speed, stability, and population diversity.
6. **Interpretation**: Draw conclusions on algorithm strengths and weaknesses across different problem types and scales.

