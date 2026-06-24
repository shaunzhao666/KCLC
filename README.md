# Introduction 
This is the Github Repository for paper *Improving Phenotype Clustering in Multivariate Phenome-wide Association Studies and the Study of Pleiotropic Effects of Susceptible Genes of Lung and Prostate Cancers*

## code
### `functions.R`: functions for simulation.

#### `lambdas_setup()`
This function generate effect size corresponding to paper descriptions
```r
lambdas_setup(beta, scene, K)
```
##### Arguments
- beta: control the effect size
- scene: six scenes referring to the paper
- K: total number of phenotypes
##### Values
effect sizes organized by category, where each category $m$ has a vector of length $k_m$:
$$\left{ \lambda_{11}, \ldots, \lambda_{1k_1} \right}, \left{ \lambda_{21}, \ldots, \lambda_{2k_2} \right}, \ldots, \left{ \lambda_{M1}, \ldots, \lambda_{Mk_M} \right}$$

####  `simulation()`
This function generates simulated phenotypic data and a genetic variant. It is designed for power analysis and method comparison in this paper.
```r
simulation(n, MAF, lambda, rho_a, rho_w, c = sqrt(0.5), seed = NA)
```
##### Arguments
- n: sample size
- MAF: minor allele frequency
- lambda: $\left\{\lambda_{11}, \ldots, \lambda_{1k_1}\right\}, \left\{\lambda_{21}, \ldots, \lambda_{2k_2}\right\}, \ldots, \left\{\lambda_{M1}, \ldots, \lambda_{Mk_M}\right\}$ effect size for each category, could be generate by `lambdas_setup(beta, scene, K)` function
- rho_a: control correlation across phenotypic categoris
- rho_w: controls the correlation among phenotypes within each category
- c: constant
- seed: optional
##### Values
Returns a list containing:
- y: phenotypes
- x: genotypes, {0, 1, 2}
- cat: information of categories based on lambda



### `script.R`: reproduce the simulation results


## simulation_results





