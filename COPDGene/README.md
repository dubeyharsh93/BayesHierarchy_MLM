# COPDGene Study Bayesian Hierarchical Model Analysis

Reproducible Julia implementation of the Bayesian hierarchical model analysis for the COPDGene Study metabolomics application presented in the accompanying manuscript.

## Table of Contents

- [Introduction](#introduction)
- [Analysis Design](#analysis-design)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Results](#results)
- [Dependencies](#dependencies)

## Introduction

This directory contains the reproducible Julia implementation of the Bayesian hierarchical model analysis for the COPDGene Study metabolomics application. The project is distributed as a self-contained Julia environment and reproduces the repeated-subsampling results and Reference MSE Ratio figure reported for COPDGene.

## Analysis Design

The analysis evaluates how the benefit of biologically informed partial pooling changes with the amount of available training data. Complementary training and testing subsets are repeatedly sampled at training sizes of `40`, `60`, `80`, `100`, `125`, `150`, `200`, `250`, `300`, `400`, `500`, and `600`. Each training size is evaluated using 100 repetitions by default.

For each split, MatrixLM is fitted independently to the training and testing subsets. The Bayesian hierarchical model is applied to the training MatrixLM effect estimates and standard errors. Because the true metabolite effects are unknown in real data, both training estimators are compared with the MatrixLM estimates from the complementary test subset.

For each covariate, the Reference MSE Ratio is

```text
mean MatrixLM reference MSE / mean Bayesian reference MSE
```

where each mean is taken across repeated splits. Values greater than one indicate closer agreement of the Bayesian training estimates with the independent reference estimates. Error Reduction reports the corresponding percentage decrease in mean reference MSE.

## Project Structure

- **BayesHierarchy/**
  - Project.toml : Julia package dependencies.
  - Manifest.toml : Pinned package versions ensuring reproducibility.
  - run_copd_mse.jl : Runs the repeated-subsampling analysis and writes covariate-specific Reference MSE summaries.
  - run_copd_plot.jl : Generates the manuscript Reference MSE Ratio figure from the saved summary table.
  - gibbs_src/ : Contains the implementation of the Bayesian hierarchical Gibbs sampler together with supporting utility functions used throughout the analysis.
  - results/ : Stores the repetition-level metrics, summary tables, and manuscript figure.

- **data/**
  - **processed/COPDGene/**: Directory for the required processed COPDGene data files.

- **notebooks/**
  - **preprocessing/**: Contains .ipynb notebooks and .jl files for data preprocessing and wrangling.

- **src/**: Source code containing functions for preprocessing and wrangling.

## Getting Started

### Prerequisites

- **Julia**: Julia 1.12.3

### Installation

Clone the repository:

```bash
git clone https://github.com/senresearch/BayesHierarchy_MLM.git
cd BayesHierarchy_MLM/COPDGene/BayesHierarchy
```

Instantiate

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

### Reproducing numerical results

```bash
julia --project=. run_copd_mse.jl
```

The default analysis uses 100 repetitions. For a one-repetition validation run:

```bash
julia --project=. run_copd_mse.jl 1
```

Produces

```text
results/tables/copdgene_rmr_all_reps.csv
results/tables/copdgene_rmr_summary.csv
results/tables/copdgene_rmr_selected_covariates.csv
```

The full 100-repetition analysis is computationally intensive. The one-repetition command should be used first to validate data paths, matrix dimensions, covariate names, and output generation.

### Reproducing manuscript figures

```bash
julia --project=. run_copd_plot.jl
```

Produces

```text
results/figures/copdgene_rmr.pdf
```

The plotting script reads the saved selected-covariate summary and does not refit either model.


## Results

All generated figures and numerical summaries are written to the `results/` directory. The manuscript figure displays Reference MSE Ratio learning curves for Age, BMI, and COPD status.

The split and Gibbs seeds are deterministic. Running the analysis with the pinned Julia environment, the same input data, and the same repetition count reproduces the same design and results.

## Dependencies

This project is distributed as a fully reproducible Julia environment. All required package versions are specified in Project.toml and Manifest.toml. Running Pkg.instantiate() automatically installs the exact package versions used in the manuscript analyses.
