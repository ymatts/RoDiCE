# RoDiCE

<!-- badges: start -->
[![R-CMD-check](https://github.com/ymatts/RoDiCE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ymatts/RoDiCE/actions/workflows/R-CMD-check.yaml)
[![test-coverage](https://github.com/ymatts/RoDiCE/actions/workflows/test-coverage.yaml/badge.svg)](https://github.com/ymatts/RoDiCE/actions/workflows/test-coverage.yaml)
[![License: GPL v3](https://img.shields.io/badge/License-GPLv3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
<!-- badges: end -->

**RoDiCE** (Robust Differential Co-Expression) is an R package for robust statistical analysis of protein complex co-expression patterns between different biological conditions (e.g., tumor vs. normal tissue).

## Overview

RoDiCE provides empirical copula-based two-sample testing for differential co-expression analysis. Unlike traditional correlation-based methods, RoDiCE is robust to noisy protein expression data and captures complex dependency structures beyond linear relationships.

### Key Features

- **Copula-based testing**: Detect differential co-expression using empirical copula statistics
- **Multivariate analysis**: Test joint dependencies across multiple proteins
- **Pairwise analysis**: Automated pairwise testing with multiple testing correction
- **Interactive visualization**: Network and scatter plot visualization of results
- **CORUM integration**: Built-in support for protein complex membership from CORUM database
- **High performance**: Parallelized C++ implementation via Rcpp

## Installation

### From GitHub (recommended)

```r
# Install devtools if not already installed
install.packages("devtools")

# Install RoDiCE
devtools::install_github("ymatts/RoDiCE")
```

### Dependencies

RoDiCE requires R >= 3.6.0 and the following packages:
- Rcpp, RcppParallel (for performance)
- gPdtest (for p-value approximation)
- ggplot2, GGally, visNetwork (for visualization)
- Rfast (for combinatorics)

## Quick Start

```r
library(RoDiCE)

# Load example data (clear cell renal cell carcinoma)
data(ccrcc.pbaf)
tumor <- ccrcc.pbaf$tumor    # 110 samples x 10 proteins
normal <- ccrcc.pbaf$normal  # 84 samples x 10 proteins

# Multivariate copula test
result <- coptest(tumor, normal, nperm = 100, approx = TRUE)
result$pval  # p-value for differential co-expression

# Pairwise copula test
result_pairwise <- coptest.p(tumor, normal, nperm = 100, approx = TRUE)
head(result_pairwise$tbl)  # Results for all protein pairs

# Visualize network
netvis(result_pairwise, title = "PBAF Complex")
```

## Main Functions

| Function | Description |
|----------|-------------|
| `coptest()` | Multivariate copula test for all joint dependencies |
| `coptest.p()` | Pairwise copula test for bivariate dependencies |
| `netvis()` | Interactive network visualization |
| `coexvis()` | Dependence pattern visualization using scatter plots |
| `group.exprs()` | Split expression data by protein complex membership |
| `parse.corum()` | Parse protein complex information from CORUM database |

## Example Workflow

```r
library(RoDiCE)

# 1. Load expression data and protein complex information
data(ccrcc.subset)
data(corum.subset)

tumor <- ccrcc.subset$tumor
normal <- ccrcc.subset$normal

# 2. Group proteins by complex
tumor_by_complex <- group.exprs(tumor, corum.subset)
normal_by_complex <- group.exprs(normal, corum.subset)

# 3. Test differential co-expression for each complex
for (i in seq_along(tumor_by_complex$comp.exprs)) {
  t_expr <- tumor_by_complex$comp.exprs[[i]]
  n_expr <- normal_by_complex$comp.exprs[[i]]

  if (ncol(t_expr) >= 2 && ncol(n_expr) >= 2) {
    result <- coptest(t_expr, n_expr, nperm = 100)
    cat(sprintf("%s: p-value = %.4f\n",
                tumor_by_complex$comp.names[i], result$pval))
  }
}
```

## Tutorial

For a detailed tutorial, please visit: https://rpubs.com/ymatts/RoDiCE

## Citation

If you use RoDiCE in your research, please cite:

> Yusuke MATSUI et al. (2020) RoDiCE: Robust differential protein co-expression analysis for cancer complexome.

## References

- Clark DJ et al. (2019) Integrated Proteogenomic Characterization of Clear Cell Renal Cell Carcinoma. Cell; 179(4), 964-983.
- Giurgiu M et al. (2019) CORUM: the comprehensive resource of mammalian protein complexes-2019. Nucleic Acids Res, 47(D1), D559-D563.

## License

GPL-3

## Author

Yusuke MATSUI (mail.to.matsui@gmail.com)
