# Simulation Bayesian Hierarchical Model Study

This directory contains the reproducible Julia implementation of the simulation study presented in the accompanying manuscript, allowing users to reproduce the performance evaluation of the Bayesian hierarchical model under controlled simulation settings.

## Table of Contents

- [Introduction](#introduction)
- [Simulation Design](#simulation-design)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Results](#results)
- [Dependencies](#dependencies)

## Introduction

The purpose of this simulation study is to evaluate the estimation and prediction performance of the proposed Bayesian hierarchical model under synthetic data generated from known hierarchical structures. The simulation reproduces the experiments presented in the manuscript and provides a controlled comparison between the proposed Bayesian approach and the corresponding Matrix Linear Model baseline.

## Simulation Design

The simulation study generates synthetic datasets from a hierarchical Bayesian model with known subclass structure. The number of metabolites is fixed at `m = 300`, while the total sample size varies over `50`, `75`, `100`, `150`, `200`, `300`, `400`, `550`, `700`, `850`, and `1000`. Each simulated dataset is divided into 70% training and 30% testing observations.

Three hierarchical heterogeneity regimes are examined:

- low: `tau_v = 0.06` and `tau_w = 0.04`;
- moderate: `tau_v = 0.12` and `tau_w = 0.08`;
- high: `tau_v = 0.24` and `tau_w = 0.16`.

The ratio `tau_v / tau_w` is held fixed at 1.5, preserving the relative allocation of between-subclass and within-subclass variability while total hierarchical heterogeneity increases. Each combination of sample size and heterogeneity regime is repeated 100 times.

Coefficient estimation performance is assessed relative to the known coefficient matrix used to generate each dataset. Prediction performance is evaluated on the independent testing observations. Both outcomes are summarized as the ratio of the mean MatrixLM MSE to the mean Bayesian MSE across repetitions, with values greater than one favoring the Bayesian estimator.

## Project Structure

- Project.toml : Julia package dependencies.
- Manifest.toml : Pinned package versions ensuring reproducibility.
- run_sim_mse.jl : Runs the complete simulation study, including data generation, Bayesian hierarchical model fitting, MatrixLM comparison, and performance evaluation.
- run_sim_plots.jl : Generates the manuscript estimation and prediction MSE ratio figures from the saved simulation summary.
- gibbs_src/ : Contains the implementation of the Bayesian hierarchical Gibbs sampler together with supporting utility functions used throughout the analysis.
- results/ : Stores all simulation outputs, including numerical summaries and intermediate results.

## Getting Started

### Prerequisites

- **Julia**: Julia 1.12.3

### Installation

Clone the repository:

```bash
git clone https://github.com/senresearch/BayesHierarchy_MLM.git
cd BayesHierarchy_MLM/Simulation
```

Instantiate

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

### Running the simulation

```bash
julia --project=. run_sim_mse.jl
```
This script

* generates synthetic datasets,
* fits both competing models,
* computes estimation and prediction errors,
* summarizes the simulation results,
* writes repetition-level and summary tables to `results/tables/`.

The full 100-repetition simulation is computationally intensive. Once it has completed, we generate the manuscript figures without rerunning the simulation:

```bash
julia --project=. run_sim_plots.jl
```

This script reads the saved three-regime summary and writes:

```text
results/figures/sim_est_ratio.pdf
results/figures/sim_pred_ratio.pdf
```

## Results

Running the simulation produces:

```text
results/tables/simulation_1level_3hetero_all_reps.csv
results/tables/simulation_1level_3hetero_summary.csv
```

The repetition-level table contains 3,300 rows corresponding to three heterogeneity regimes, eleven sample sizes, and 100 repetitions. The summary table contains the mean estimation and prediction MSEs for both methods and their corresponding ratios of mean MSEs.

## Dependencies

This project is distributed as a fully reproducible Julia environment. All required package versions are specified in Project.toml and Manifest.toml. Running Pkg.instantiate() automatically installs the exact package versions used in the manuscript analyses.
