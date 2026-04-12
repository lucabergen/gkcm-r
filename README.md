# gkcm

`gkcm` implements the [**Generalised Kernel Covariance Measure (GKCM)**](https://arxiv.org/abs/2604.03721) for conditional independence testing in R.

The package provides a regression-based, kernel-based conditional independence test that supports **mixed data**, including continuous and categorical variables. It also includes an interface for use in **constraint-based causal discovery** workflows.

## Features

- Conditional and marginal independence testing
- Support for **mixed data**:
  - Gaussian kernel for numeric variables
  - Dirac kernel for categorical variables
- Regression backend based on:
  - `drf` for numeric targets
  - `ranger` for categorical targets
- Integration with causal discovery workflows via:
  - `gkcm_suffStat()`
  - `gkcm_indepTest()`
- Optional parallel computation via `num_threads`

## Installation

You can install the development version from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("lucabergen/gkcm-r")
