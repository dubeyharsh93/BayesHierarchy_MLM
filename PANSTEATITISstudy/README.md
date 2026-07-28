# PANSTEATITIS Study Bayesian Hierarchical Model Analysis

Reproducible Julia implementation of the Bayesian hierarchical model analysis for the PANSTEATITIS metabolomics application presented in the accompanying manuscript.

## Table of Contents

- [Introduction](#introduction)
- [Analysis Design](#analysis-design)
- [Project Structure](#project-structure)
- [Required Data](#required-data)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Results](#results)
- [Dependencies](#dependencies)

## Introduction

This directory contains the reproducible Julia implementation of the Bayesian hierarchical model analysis for the PANSTEATITIS metabolomics application presented in the accompanying manuscript. The project is distributed as a self-contained Julia environment and reproduces the repeated-subsampling table reported for the Weight covariate.

## Analysis Design

The PANSTEATITIS analysis evaluates the benefit of biologically informed partial pooling in metabolite effect estimation across training sizes of `10`, `15`, `20`, and `25` of the PANSTEATITIS data. Each training size is evaluated using 100 repeated complementary training/testing splits by default.

For each split, MatrixLM is fitted independently to the training and testing subsets. The Bayesian hierarchical model is applied to the training MatrixLM effect estimates and standard errors. Both training estimators are compared with the MatrixLM effect estimates from the complementary reference subset.

For each covariate, the Reference MSE Ratio is

```text
mean MatrixLM reference MSE / mean Bayesian reference MSE
```

where each mean is taken across repeated splits. Values greater than one indicate closer agreement of the Bayesian training estimates with the independent reference estimates. Error Reduction reports the corresponding percentage decrease in mean reference MSE.

The manuscript table reports results for the `Weight` coefficient. The analysis excludes the intercept and uses curated biochemical subclass annotations for the metabolite hierarchy.

## Project Structure

- **BayesHierarchy/**
  - Project.toml : Julia package dependencies.
  - Manifest.toml : Pinned package versions ensuring reproducibility.
  - run_pans_mse.jl : Runs the repeated-subsampling analysis, writes the Reference MSE tables, and prints the manuscript table.
  - gibbs_src/ : Contains the implementation of the Bayesian hierarchical Gibbs sampler together with supporting utility functions used throughout the analysis.
  - results/ : Stores repetition-level metrics and summary tables.

- **data/**
  - **processed/**: Directory for the required processed PANSTEATITIS data files.

- **notebooks/**
  - **preprocessing/**: Contains notebooks and Julia files for data preprocessing and wrangling.

- **src/**: Source code containing functions for preprocessing and wrangling.

## Getting Started

### Prerequisites

- **Julia**: Julia 1.12.3

### Installation

Clone the repository:

```bash
git clone https://github.com/senresearch/BayesHierarchy_MLM.git
cd BayesHierarchy_MLM/PANSTEATITISstudy/BayesHierarchy
```

Instantiate the pinned Julia environment:

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

### Reproducing numerical results

Run the complete analysis:

```bash
julia --project=. run_pans_mse.jl
```

The default analysis uses 100 repetitions. For a one-repetition validation run:

```bash
julia --project=. run_pans_mse.jl 1
```

The analysis produces:

```text
results/tables/pansteatitis_rmr_all_reps.csv
results/tables/pansteatitis_rmr_summary.csv
results/tables/pansteatitis_rmr_weight.csv
```

The script also prints a booktabs-formatted LaTeX version of the manuscript table directly in the terminal. The full 100-repetition analysis is computationally intensive, so the one-repetition command should be used first to validate data paths, matrix dimensions, covariate names, and output generation.

## Results

All generated numerical summaries are written to `results/tables/`.

The split and Gibbs seeds are deterministic. Running the analysis with the pinned Julia environment, the same input data, and the same repetition count reproduces the same design and results.

## Dependencies

This project is distributed as a fully reproducible Julia environment. All required package versions are specified in `Project.toml` and `Manifest.toml`. Running `Pkg.instantiate()` automatically installs the exact package versions used in the manuscript analysis.
