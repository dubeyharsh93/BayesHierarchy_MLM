# A Bayesian framework using shared features improves metabolite estimation and prediction

Harsh Vardhan Dubey<sup>1</sup>, Gregory Farage<sup>1</sup>, Katerina Kechris<sup>2</sup>, Śaunak Sen<sup>1</sup>

><sup>1</sup>Department of Preventive Medicine, College of Medicine, University of Tennessee Health Science Center, Memphis, TN   
<sup>2</sup>Department of Biostatistics & Informatics, Colorado School of Public Health, University of Colorado Anschutz Medical Campus, Aurora, CO

### Abstract

In metabolomics studies, incorporating shared biochemical feature information 
among metabolites can improve estimation and prediction of metabolite effects in high-dimensional settings with limited sample sizes. 
Many existing approaches treat metabolites as independent features, ignoring known biochemical structure such as shared subclasses and pathway membership. 
We introduce a Bayesian hierarchical model that improves effect estimates from any well-established analytical method, using MatrixLM as a baseline in our evaluation, by incorporating the multilevel organization of 
metabolites into subclasses and broader biochemical categories. 
Bayesian shrinkage stabilizes individual metabolite estimates through partial pooling, reducing mean squared error while preserving interpretability. 
We evaluate the approach using simulation studies across varying sample sizes and heterogeneity regimes, along with applications to three metabolomics datasets. 
Our analyses suggest that incorporating shared feature 
information improves both estimation and prediction relative to non-Bayesian linear 
alternatives.

### Repository Info

The repository reproduces all analyses presented in the manuscript, including

* three real metabolomics applications
* one simulation study
* manuscript figures
* manuscript tables

Each application is distributed as an independent Julia project with its own Project.toml and Manifest.toml, ensuring fully reproducible computational environments.

### Requirements

* Julia 1.12.3

Each project includes a fully specified Julia environment through

* Project.toml
* Manifest.toml

No additional package installation is required.

### Reproducing the Manuscript

Clone the repository

```bash
git clone https://github.com/<username>/BayesHierarchy_MLM.git
cd BayesHierarchy_MLM
```

Navigate to the project you wish to reproduce.

For example,
```bash
cd COPDGene/BayesHierarchy
```

Instantiate the Julia environment
```bash
julia --project=. -e 'using Pkg; Pkg.instantiate()'
```

This installs all package versions required for the analysis.

#### COPDGene
```bash
cd COPDGene/BayesHierarchy
julia run_copd_mse.jl
julia run_copd_plot.jl
```

#### SAMS
```bash
cd SAMSstudy/BayesHierarchy
julia run_sams_mse.jl
julia run_sams_plot.jl
```

#### PANSTEATITIS
```bash
cd PANSTEATITISstudy/BayesHierarchy
julia run_pans_mse.jl
julia run_pans_plot.jl
```

#### Simulation

```bash
cd Simulation/BayesHierarchy
julia run_sim_mse.jl
```

### Expected Outputs

Running the scripts reproduces

* manuscript figures
* manuscript summary tables
* posterior estimates
* prediction accuracy measures
* estimation reproducibility analyses

All generated outputs are written to the corresponding results/ directory.

### Citation

If you use this repository, please cite

> Paper citation here

