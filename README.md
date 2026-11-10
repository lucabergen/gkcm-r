# gkcm

`gkcm` implements the **Generalised Kernel Covariance Measure (GKCM)** for conditional independence testing in R.

The package provides a kernel-based conditional independence test that supports **mixed data**, including continuous and categorical variables. It also includes an interface for use in **constraint-based causal discovery** workflows.

## Features

`gkcm` supports both conditional and marginal independence testing. By default, random Fourier features approximate the Gaussian-kernel feature map for
numeric variables, helping GKCM scale to larger samples.
Forest regression runs in parallel; with the default
`num_threads = NULL`, half of the available CPU cores are used. Increase `num_threads` to
use more cores and reduce runtime. RKHS-valued regression uses `okranger` and `ranger`;
`gkcm_suffStat()` and `gkcm_indepTest()` provide interfaces for constraint-based causal
discovery.

## Installation

Install from GitHub with:

```r
# install.packages("remotes")
remotes::install_github("lucabergen/gkcm-r")
```

## Quick start

```r
library(gkcm)
set.seed(1)

n <- 200
z <- data.frame(z = rnorm(n))
x <- z$z + rnorm(n, sd = 0.5)
y <- z$z + rnorm(n, sd = 0.5)

# Conditional test
gkcm(x, y, z = z)

# Marginal test
gkcm(x, y, z = NULL)

# Categorical outcomes must be factors
y_cat <- factor(ifelse(x > median(x), "high", "low"))
gkcm(x, y_cat, z = z)
```

## Citation

Please cite the GKCM paper when using this package:

- Bergen, L., Sejdinovic, D., & Didelez, V. (2026). [*The Generalised Kernel Covariance Measure*](https://arxiv.org/abs/2604.03721). *Proceedings of Machine Learning Research*, 323, 1–23.

Related references on random Fourier features and their use in conditional independence testing:

- Rahimi, A., & Recht, B. (2007). [*Random Features for Large-Scale Kernel Machines*](https://papers.nips.cc/paper_files/paper/2007/hash/013a006f03dbc5392effeb8f18fda755-Abstract.html). *Advances in Neural Information Processing Systems*, 20.
- Strobl, E. V., Zhang, K., & Visweswaran, S. (2019). [*Approximate Kernel-Based Conditional Independence Tests for Fast Non-Parametric Causal Discovery*](https://doi.org/10.1515/jci-2018-0017). *Journal of Causal Inference*, 7(1).
