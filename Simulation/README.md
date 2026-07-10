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

The simulation study generates synthetic datasets from a hierarchical Bayesian model with known subgroup structure. The generated datasets are used to evaluate the estimation and prediction performance of the proposed Bayesian hierarchical model relative to the Matrix Linear Model baseline under repeated Monte Carlo replications.

## Project Structure

- Project.toml : Julia package dependencies.
- Manifest.toml : Pinned package versions ensuring reproducibility.
- run_sim_mse.jl : Runs the complete simulation study, including data generation, Bayesian hierarchical model fitting, MatrixLM comparison, and performance evaluation.
- gibbs_src/ : Contains the implementation of the Bayesian hierarchical Gibbs sampler together with supporting utility functions used throughout the analysis.
- results/ : Stores all simulation outputs, including numerical summaries and intermediate results.

## Getting Started

### Prerequisites

- **Julia**: Julia 1.12.3

### Installation

Clone the repository:

```bash
git clone https://github.com/<yourname>/BayesHierarchy_MLM.git
cd BayesHierarchy_MLM/Simulation
```

Instantiate

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

### Running the simulation

```bash
julia run_sim_mse.jl
```
This script

* generates synthetic datasets,
* fits both competing models,
* computes estimation and prediction errors,
* summarizes the simulation results,
* writes the outputs to the results/ directory.

## Results

Running the simulation reproduces the numerical results reported in the manuscript. Since the study relies on Monte Carlo sampling, small numerical differences (typically on the order of 10^{-3}) relative to the values in the paper are expected and do not affect the overall conclusions.

## Dependencies

This project is distributed as a fully reproducible Julia environment. All required package versions are specified in Project.toml and Manifest.toml. Running Pkg.instantiate() automatically installs the exact package versions used in the manuscript analyses.
