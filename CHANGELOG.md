# Changelog

All notable changes to RoDiCE will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added
- **Testing framework**: Added testthat framework with comprehensive test coverage
  - Tests for `coptest()` function
  - Tests for `coptest.p()` function
  - Tests for `group.exprs()` function
  - Tests for utility functions
- **CI/CD pipeline**: Added GitHub Actions workflows
  - R CMD check on multiple platforms (Linux, macOS, Windows)
  - Test coverage reporting
- **Utility functions**: New `R/utils.R` file with reusable functions
  - `gpd_pval_approx()`: GPD-based p-value approximation
  - `validate_matrix_input()`: Input validation for matrix arguments
- **Documentation**: Comprehensive README.md with
  - Package overview and features
  - Installation instructions
  - Quick start guide
  - Example workflow
  - Function reference table

### Changed
- **Code refactoring**: Extracted duplicated GPD approximation logic into shared function
- **netvis.R**: Refactored into modular helper functions
  - `.validate_netvis_input()`: Input validation
  - `.style_nodes()`: Node styling
  - `.style_links()`: Link styling
  - `.create_legend()`: Legend creation
  - `.netvis_without_ref()`: Visualization without reference
  - `.netvis_with_ref()`: Visualization with reference
- **coptest.R**: Added input validation and simplified GPD logic
- **coptest_p.R**:
  - Replaced `print()` with `message()` for progress output
  - Added proper input validation

### Fixed
- **Bug fix in coptest_p.R**: Fixed incorrect null distribution reference in GPD approximation loop (was using stale `null` variable instead of `result_perm[[i]]$null`)
- **Input validation**: Added comprehensive validation for matrix inputs
  - NULL checks
  - Type checks (matrix, numeric)
  - Dimension checks
  - NA/NaN/Inf checks
  - Column name validation

### Improved
- **Error messages**: More descriptive and actionable error messages throughout
- **Code style**: Consistent formatting and coding style
- **Documentation**: Added roxygen2 documentation for internal functions

## [0.99.1] - Previous Release

### Features
- Initial implementation of copula-based differential co-expression testing
- `coptest()`: Multivariate copula test
- `coptest.p()`: Pairwise copula test
- `netvis()`: Network visualization
- `coexvis()`: Co-expression visualization
- `group.exprs()`: Expression data grouping
- `parse.corum()`: CORUM database parsing
- Parallelized C++ implementation via Rcpp/RcppParallel
- Example datasets from CCRCC study
