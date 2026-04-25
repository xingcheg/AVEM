# AVEM (Anchored Variational EM)

This repository collects the corresponding R files and saved outputs for AVEM-based simulation studies in two model families:

- `AVEM_HMM`: mean-anchor AVEM for Gaussian mixed hidden Markov models (MHMMs).
- `AVEM_SSM`: AVEM for mixed-effects state-space models (MESSM).

The code is organized as standalone R scripts rather than as an R package, so the main workflow is to `source()` the scripts and run the experiment drivers directly.

## Repository Structure

```text
.
|-- AVEM_HMM/
|   |-- main.R
|   |-- experiment1.R
|   |-- experiment2.R
|   |-- gaussian_mhmm_compare_methods.R
|   |-- bernoulli_mhmm_compare_methods.R
|   |-- accuracy_experiment.rds
|   `-- speed_experiment.rds
`-- AVEM_SSM/
    |-- avem_messm_simulation.R
    |-- avem_messm_mc.R
    |-- experiment1.R
    |-- visualize_exp1.R
    |-- visualize_exp1_elbo.R
    `-- results/
        |-- case1_out.rds
        |-- case2_out.rds
        |-- case3_out.rds
        |-- case4_out.rds
        |-- case5_out.rds
        `-- case6_out.rds
```

## Directory Guide

### `AVEM_HMM`

This folder contains code for Gaussian mixed hidden Markov models with a common latent Markov chain and subject-level random effects.

- `main.R`: core implementation for Gaussian MHMM simulation and mean-anchor AVEM fitting.
- `experiment1.R`: accuracy study that varies sample size `n`, sequence length `T_len`, and random-effect variance `tau2`.
- `experiment2.R`: speed study that varies the number of latent states `K` and observation dimension `d`.
- `gaussian_mhmm_compare_methods.R`: comparison script for Gaussian MHMMs, including mean-anchor AVEM, Gaussian quadrature EM, and Monte Carlo EM.
- `bernoulli_mhmm_compare_methods.R`: comparison script for Bernoulli/logistic MHMMs, including Laplace-AVEM, Gaussian quadrature EM, and Monte Carlo EM.
- `accuracy_experiment.rds`: saved output from the Gaussian MHMM accuracy experiment.
- `speed_experiment.rds`: saved output from the Gaussian MHMM speed experiment.

### `AVEM_SSM`

This folder contains code for mixed-effects state-space models with subject-specific transition and loading parameters.

- `avem_messm_simulation.R`: self-contained implementation of data generation, Kalman smoothing, AVEM updates, and evaluation metrics.
- `avem_messm_mc.R`: Monte Carlo wrapper that stores replicate-level metrics and ELBO trajectories.
- `experiment1.R`: example Monte Carlo experiment using the SSM implementation.
- `visualize_exp1.R`: summary boxplots for RMSE-based metrics across saved experiment cases.
- `visualize_exp1_elbo.R`: ELBO trajectory visualization across saved experiment cases.
- `results/case*_out.rds`: saved Monte Carlo outputs used by the visualization scripts.

## Requirements

The scripts use base R plus a small set of plotting and data-manipulation packages.

Install the required packages with:

```r
install.packages(c("ggplot2", "gridExtra", "dplyr", "tidyr", "purrr", "statmod"))
```

## Quick Start

### 1. HMM example

```r
setwd("AVEM_HMM")
source("main.R")

res <- run_one_replication(
  K = 3,
  n = 50,
  T_len = 60,
  d = 2,
  tau2 = 1,
  seed = 2026
)

print(res)
```

To run the bundled studies:

```r
source("experiment1.R")  # accuracy experiment
source("experiment2.R")  # speed experiment
```


#### Comparing AVEM with QEM and MCEM

For Gaussian mixed HMMs:

```r
setwd("AVEM_HMM")
source("gaussian_mhmm_compare_methods.R")

one_gaussian <- run_one_replication_compare(
  K = 2,
  n = 40,
  T_len = 40,
  d = 2,
  tau2 = 1,
  gh_nodes = c(3, 5),
  mc_points = c(25, 50),
  seed = 2026,
  verbose = TRUE
)

print(one_gaussian$results)
```

For Bernoulli/logistic mixed HMMs:

```{r}
setwd("AVEM_HMM")
source("bernoulli_mhmm_compare_methods.R")

one_bernoulli <- run_one_replication_compare_bernoulli(
  n = 40,
  T_len = 100,
  beta_true = c(-1.5, 1.5),
  tau2 = 1,
  gh_nodes = c(10, 20),
  mc_points = c(10, 20),
  seed = 2026,
  verbose = TRUE
)

print(one_bernoulli$results)
```

### 2. SSM example

```r
setwd("AVEM_SSM")
source("avem_messm_simulation.R")
source("avem_messm_mc.R")

out <- run_one_simulation(config)
```

For a Monte Carlo run:

```r
base_config <- config
base_config$seed <- 2026
base_config$n <- 50
base_config$T <- 100
base_config$p <- 4
base_config$q <- 2
base_config$max_iter <- 120
base_config$rel_tol <- 1e-4

mc_out <- run_monte_carlo_avem_all_curves(
  base_config = base_config,
  n_rep = 100,
  seed_start = 2026
)
```

To reproduce the supplied plots:

```r
source("visualize_exp1.R")
source("visualize_exp1_elbo.R")
```

## Saved Results

This repository includes precomputed `.rds` files so that figures and summaries can be reproduced without rerunning every simulation:

- `AVEM_HMM/accuracy_experiment.rds`
- `AVEM_HMM/speed_experiment.rds`
- `AVEM_SSM/results/case1_out.rds` through `case6_out.rds`

These files are especially useful because some of the example studies use `n_rep = 100` and may take noticeable time to rerun.

