# BayesHierarchy_MLM
Bayesian Hierarchical Model for Metabolomic Studies

# Model Overview

## Bilinear Regression Model (MatrixLM)

Let  
$Y \in \mathbb{R}^{n \times m}$ denote a matrix of responses (e.g., metabolite abundances),  
$X \in \mathbb{R}^{n \times p}$ a design matrix of covariates, and  
$Z \in \mathbb{R}^{m \times q}$ a feature design matrix encoding structure among metabolites  
(e.g., identity, subclasses, or superclasses).

The bilinear model is

$$
Y = X B Z^\top + E,
$$

where

- $B \in \mathbb{R}^{p \times q}$ is the coefficient matrix,
- $E$ is a noise matrix with independent rows.

MatrixLM provides, for each covariate $k$ and metabolite $j$,

$$
\hat{b}_{kj}, \qquad \mathrm{se}_{kj}.
$$

For a fixed covariate $k$, we treat

$$
\hat{b}_j = \hat{b}_{kj}, \qquad
s_j = \mathrm{se}_{kj}, \qquad
j = 1,\dots,m,
$$

as noisy observations of latent true effects $\theta_j$:

$$
\hat{b}_j \mid \theta_j \sim \mathcal{N}(\theta_j, s_j^2).
$$

---

## Bayesian Hierarchical Meta-Analysis Model

Metabolites are organized into subclasses and superclasses:

$$
j \;\mapsto\; h(j) \in \{1,\dots,H\}, \qquad
h \;\mapsto\; g(h) \in \{1,\dots,G\}.
$$

The hierarchy is

$$
\begin{aligned}
\hat{b}_j \mid \theta_j &\sim \mathcal{N}(\theta_j, s_j^2), \\[4pt]
\theta_j \mid \beta_{h(j)}, \tau_w^2 &\sim \mathcal{N}(\beta_{h(j)}, \tau_w^2), \\[4pt]
\beta_h \mid \alpha_{g(h)}, \tau_v^2 &\sim \mathcal{N}(\alpha_{g(h)}, \tau_v^2), \\[4pt]
\alpha_g \mid \theta_0, \tau_u^2 &\sim \mathcal{N}(\theta_0, \tau_u^2), \\[6pt]
\theta_0 &\sim \mathcal{N}(\mu_0, s_0^2).
\end{aligned}
$$

Shrinkage priors on the scale parameters:

$$
\tau_w \sim \text{Half-Cauchy}(1), \quad
\tau_v \sim \text{Half-Cauchy}(1), \quad
\tau_u \sim \text{Half-Cauchy}(1).
$$

Each Half-Cauchy is implemented via an inverse-gamma scale mixture:

$$
\tau^2 \mid \lambda \sim \text{IG}\!\left(\tfrac{1}{2}, \tfrac{1}{\lambda}\right),
\qquad
\lambda \sim \text{IG}\!\left(\tfrac{1}{2}, 1\right),
$$

yielding conditionally conjugate updates.