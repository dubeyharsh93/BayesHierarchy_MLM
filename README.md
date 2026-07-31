# A Bayesian framework using shared features improves metabolite estimation and prediction

Harsh Vardhan Dubey<sup>1</sup>, Gregory Farage<sup>1</sup>, Śaunak Sen<sup>1</sup>

><sup>1</sup>Department of Preventive Medicine, College of Medicine, University of Tennessee Health Science Center, Memphis, TN   

### Abstract

External biological knowledge provides valuable information about relationships among metabolites, yet this information is usually not incorporated directly into statistical estimation procedures.
Most existing approaches estimate metabolite effects independently, ignoring known biochemical structure such as shared subclasses and pathway membership.
We propose a Bayesian hierarchical framework that improves metabolite effect estimates by incorporating external biological information describing relationships among metabolites.
The proposed method improves metabolite-specific estimates by allowing related metabolites to borrow information from one another while preserving metabolite-level inference.
We evaluate the methodology using simulation studies across a range of sample sizes and heterogeneity regimes together with three metabolomics applications involving distinct biological annotation structures.
Across both simulated and real datasets, incorporating external biological information consistently improves metabolite effect estimation.
Gains are most pronounced when sample sizes are small and metabolite classes are informative, i.e. more homogenous within classes.

### Materials

The repository reproduces all analyses presented in the manuscript, including

* three real metabolomics applications
* one simulation study
* manuscript figures
* manuscript tables

Each application is distributed as an independent Julia project with its own Project.toml and Manifest.toml, ensuring fully reproducible computational environments.
- [**COPDGene Study**](https://github.com/senresearch/BayesHierarchy_MLM/tree/main/COPDGene)
- [**Statin-Associated Muscle Symptoms (SAMS) Study**](https://github.com/senresearch/BayesHierarchy_MLM/tree/main/PANSTEATITISstudy)  
- [**Pansteatitis Mozambique Tilapia Study**](https://github.com/senresearch/BayesHierarchy_MLM/tree/main/SAMSstudy)
- [**Simulation**](https://github.com/senresearch/BayesHierarchy_MLM/tree/main/Simulation) 

### Data availability

Study data are not distributed with this repository. Each application README lists the required input files and their expected locations. After obtaining the relevant data, place the files in the specified data directory before running the corresponding reproducibility pipeline.

### Citation

If you use this repository, please cite

- Dubey, H. V., et al. Using Shared Features Improves Metabolite Effect Estimation. bioRxiv preprint, 2026.

### References:

- Gregory Farage, Chenhao Zhao, Hyo Young Choi, Timothy J. Garrett, Marshall B. Elam, Katerina Kechris, and Śaunak Sen. Matrix linear models for connecting metabolite composition to individual characteristics. Metabolites, 15(2), 2025. ISSN 2218-1989. doi: 10.3390/metabo15020140. URL https://www.mdpi.com/2218-1989/15/2/140.

### Resources:

- [MatrixLM.jl package](https://github.com/senresearch/MatrixLM.jl)
- [Sen Research Group Resources](https://senresearch.github.io/)