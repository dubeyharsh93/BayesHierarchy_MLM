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

### Materials

The repository reproduces all analyses presented in the manuscript, including

* three real metabolomics applications
* one simulation study
* manuscript figures
* manuscript tables

Each application is distributed as an independent Julia project with its own Project.toml and Manifest.toml, ensuring fully reproducible computational environments.
- [**COPDGene Study**](https://github.com/dubeyharsh93/BayesHierarchy_MLM/tree/main/COPDGene)
- [**Statin-Associated Muscle Symptoms (SAMS) Study**](https://github.com/dubeyharsh93/BayesHierarchy_MLM/tree/main/PANSTEATITISstudy)  
- [**Pansteatitis Mozambique Tilapia Study**](https://github.com/dubeyharsh93/BayesHierarchy_MLM/tree/main/SAMSstudy)
- [**Simulation**](https://github.com/dubeyharsh93/BayesHierarchy_MLM/tree/main/Simulation)  

### Citation

If you use this repository, please cite

> Paper citation here

### References:

- Gregory Farage, Chenhao Zhao, Hyo Young Choi, Timothy J. Garrett, Marshall B. Elam, Katerina Kechris, and Śaunak Sen. Matrix linear models for connecting metabolite composition to individual characteristics. Metabolites, 15(2), 2025. ISSN 2218-1989. doi: 10.3390/metabo15020140. URL https://www.mdpi.com/2218-1989/15/2/140.

### Resources:

- [MatrixLM.jl package](https://github.com/senresearch/MatrixLM.jl)
- [Sen Research Group Resources](https://senresearch.github.io/)