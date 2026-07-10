# SAMS Study Bayesian Hierarchical Model Analysis

Reproducible Julia implementation of the Bayesian hierarchical model analysis for the SAMS Study metabolomics application presented in the accompanying manuscript.

## Table of Contents

- [Introduction](#introduction)
- [Project Structure](#project-structure)
- [Getting Started](#getting-started)
- [Usage](#usage)
- [Results](#results)
- [Dependencies](#dependencies)

## Introduction

This directory contains the reproducible Julia implementation of the Bayesian hierarchical model analysis for the SAMS study metabolomics application presented in the accompanying manuscript. The project is distributed as a self-contained Julia environment and reproduces the manuscript figures, posterior estimates, and numerical summaries for the SAMS study.

## Project Structure

- **BayesHierarchy/**
  - Project.toml : Julia package dependencies.
  - Manifest.toml : Pinned package versions ensuring reproducibility.
  - run_sams_mse.jl : Runs the Bayesian hierarchical model on the SAMS metabolomics data and reproduces the estimation and prediction performance reported in the manuscript.
  - run_sams_plot.jl : Generates the manuscript figures for the SAMS application.
  - gibbs_src/ : Contains the implementation of the Bayesian hierarchical Gibbs sampler together with supporting utility functions used throughout the analysis.
  - results/ : Stores all generated output including posterior summaries, prediction results, figures, and intermediate outputs.

- **data/**
  - **processed/**: Directory for processed data files.

- **notebooks/**
  - **preprocessing/**: Contains .ipynb notebooks and .jl files for data preprocessing and wrangling.

- **src/**: Source code containing functions for preprocessing and wrangling.

## Getting Started

### Prerequisites

- **Julia**: Julia 1.12.3

### Installation

Clone the repository:

```bash
git clone https://github.com/<yourname>/BayesHierarchy_MLM.git
cd BayesHierarchy_MLM/SAMSstudy/BayesHierarchy
```

Instantiate

```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

## Usage

### Reproducing numerical results

```bash
julia run_sams_mse.jl
```
Produces

* posterior estimates
* prediction summaries
* estimation reproducibility metrics

### Reproducing manuscript figures

```bash
julia run_sams_plot.jl
```

Produces

* SAMS manuscript figures


## Results

All generated figures, numerical summaries, and intermediate outputs are written to the results/ directory.

## Dependencies

This project is distributed as a fully reproducible Julia environment. All required package versions are specified in Project.toml and Manifest.toml. Running Pkg.instantiate() automatically installs the exact package versions used in the manuscript analyses.
