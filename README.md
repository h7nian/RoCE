# FACE-HD: Federated Adaptive Causal Estimation in High Dimensions

[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Naming convention**: The paper refers to the method as **FACE-HD**. The R package
> is named `FACEHD` (R packages cannot contain hyphens).

## Overview

FACE-HD implements federated learning algorithms for causal inference in multi-site settings. The package provides communication-efficient protocols with two-level cross-fitting for enhanced robustness against overfitting bias.

### Target Estimand

**Important**: This package estimates the **potential outcome mean** at the target site:

$$\mu^1_t = E_t[Y(1)]$$

This is the expected outcome under treatment (A=1) for the **target site population**. To estimate the Average Treatment Effect (ATE), run the algorithm twice:

$$\text{ATE} = \mu^1_t - \mu^0_t = E_t[Y(1)] - E_t[Y(0)]$$

The framework extends symmetrically to estimate μ⁰ = E_t[Y(0)] by setting `A_val = 0`.

### Key Features

- **Doubly Robust Estimators**: Calibrated loss functions for Neyman Orthogonality
- **Two Communication Protocols**: Two-round and one-round algorithms
- **Two-Level Cross-fitting**: Enhanced robustness through nested sample splitting
- **Optimal Aggregation**: Variance-minimizing combination of site-specific estimators
- **C++ Acceleration**: High-performance coordinate descent with GLMNET-style optimizations
- **Comparison Methods**: Target-only, sample-size weighted, inverse-variance weighted, and tilted AIPW estimators

### API Stability and Compatibility

- **Canonical cross-fitting API**: `run_crossfit()` with `communication_mode = "two_round"` or `"one_round"`.
- **Internal C++ bindings**: functions ending with `_cpp` are internal implementation/testing interfaces and are **not** part of the package's stable public API contract.

## Installation

### Prerequisites

Ensure you have the following R packages installed:

```r
install.packages(c("Rcpp", "RcppEigen", "glmnet", "MASS", "parallel", "doParallel", "foreach"))

# For comparison methods (optional)
install.packages("RCAL")

# For testing (optional)
install.packages("testthat")
```

### From Source

```bash
git clone https://github.com/sinianzhang/FACE-HD.git
cd FACE-HD
```

In R:

```r
# Option 1: Install the package
# R CMD INSTALL .

# Option 2: Load for development (requires devtools)
devtools::load_all(".")
```

## Quick Start

```r
# Generate simulation data
set.seed(42)
data <- generate_simulation_data(
  n_total = 5000,  # Total sample size
  K = 3,           # Number of source sites
  p = 10,          # Number of covariates
  config = "C1"    # Configuration (C1-C4)
)

# Split data by site
data_split <- split_data_by_site(data)

# Run two-round cross-fitting algorithm (canonical API)
# NOTE: This estimates μ¹_t = E_t[Y(1)], NOT ATE
result_mu1 <- run_crossfit(
  data_split,
  communication_mode = "two_round",
  n_folds = 10,             # Cross-fitting folds
  lambda_selection = "cv",  # Lambda selection method
  verbose = TRUE,
  M_tau = 10.0,             # Truncation parameter
  n_cores = -1              # Parallel: -1 = all cores minus 1, NULL = sequential
)

# Access results for μ¹_t = E_t[Y(1)]
cat("μ¹ Estimate:", result_mu1$estimate, "\n")
cat("SE:", result_mu1$se, "\n")
cat("True μ¹:", data$mu1_true, "\n")

# To compute ATE, also estimate μ⁰_t = E_t[Y(0)]
# (requires running with A_val = 0 in fit_initial_outcome calls)
# ATE = result_mu1$estimate - result_mu0$estimate
```

## Project Structure

```
FACE-HD/
├── R/                              # Core R functions
│   ├── constants.R                 # Centralized numerical constants
│   ├── model_fitting.R             # Unified model fitting (R + C++ accelerated)
│   ├── cross_fitting_algorithms.R  # Two-round and one-round cross-fitting algorithms
│   ├── cross_fitting_aggregation.R # Shared aggregation utilities (DRY)
│   ├── comparison_methods.R        # Baseline methods: target-only, SS, IVW, Tilted-AIPW
│   ├── data_generation.R           # Simulation data generation (C1-C4 configs)
│   ├── data_generation_face.R      # FACE paper DGP
│   ├── estimators_oracle.R         # Oracle DR estimator (uses true parameters)
│   ├── simulation.R                # Monte Carlo simulation driver
│   ├── numerical_utils.R           # Pure numerical/math functions
│   ├── validation.R                # Input validation functions
│   └── checkpoint.R                # Checkpoint/restart for SLURM preemption
├── src/                            # C++ source files (Rcpp/RcppEigen)
│   ├── density_ratio.cpp           # Density ratio fitting and CV
│   ├── outcome_model.cpp           # Outcome model fitting, CV, and GLM utilities
│   ├── variance.cpp                # Variance, covariance, and correction terms
│   ├── weight_optimization.cpp     # Weight optimization and aggregation
│   ├── optimization.hpp            # C++ header declarations
│   ├── cv_utils.hpp                # Cross-validation shared utilities
│   ├── numerical_constants.hpp     # Numerical constants (C++ side)
│   └── utils.hpp                   # GLM utilities and numerical helpers
├── tests/                          # Test suite (testthat)
│   ├── testthat.R                  # Test runner
│   └── testthat/                   # Individual test files
├── main.R                          # Main simulation script (entry point)
├── main.sh                         # SLURM batch job submission
├── main.cmd                        # SLURM job configuration with preemption handling
├── docs/
│   ├── main.tex                    # Mathematical methodology document
│   └── proof.tex                   # Appendix proofs and assumption mapping
├── DESCRIPTION                     # R package description
├── LICENSE                         # MIT License
└── README.md                       # This file
```

## Algorithms

### Two-Round Communication Protocol

Following the methodology in `docs/main.tex` (Section A.2), the two-round algorithm:

1. **Round 1**: Source sites compute initial outcome models and send to target
2. **Target Response**: Target computes summary statistics and sends back
3. **Round 2**: Sources compute calibrated parameters and final estimates

### One-Round Communication Protocol

A more communication-efficient variant:

1. **Single Round**: Target computes all initial models and summaries
2. **Sources Respond**: Each source computes calibrated parameters and estimates

### Two-Level Cross-fitting

Both algorithms use nested sample splitting (K_f folds) to achieve:
- Neyman Orthogonality through calibrated loss functions
- Cross-calibrated plug-in order where calibrated α uses γ_init as the weight plug-in
- Unbiased estimation even with data-adaptive nuisance estimation
- Theoretical variance formula: $\text{Var}(\hat{\mu}_{agg}) = \frac{1}{K_f^2} \sum_{k_1=1}^{K_f} \text{Var}(\hat{\mu}_{agg,k_1})$

## Simulation Configurations

The FACE-HD DGP uses $X^\dagger$ for the true site/treatment and outcome
mechanisms. Configurations change only the fitted working bases exposed to the
estimators:

| Config | Fitted Site Basis | Fitted Outcome Basis | Description |
|--------|-------------------|----------------------|-------------|
| C1 | $X^\dagger$ | $X^\dagger$ | Both correctly specified |
| C2 | $X^\dagger$ | $X$ | Outcome model misspecified |
| C3 | $X$ | $X^\dagger$ | Propensity score misspecified |
| C4 | $X$ | $X$ | Both misspecified |

## Running Simulations on HPC (SLURM)

```bash
# Submit jobs for multiple parameter combinations
./main.sh
```

The simulation supports checkpointing for long-running jobs with SLURM preemption handling.

## Testing

```r
# Run all tests
testthat::test_dir("tests/testthat")

# Run specific test file
testthat::test_file("tests/testthat/test-utils.R")
```

## Citation

If you use this package in your research, please cite:

```bibtex
@software{facehd2025,
  author = {Zhang, Sinian},
  title = {FACE-HD: Federated Adaptive Causal Estimation in High Dimensions},
  year = {2025},
  url = {https://github.com/sinianzhang/FACE-HD}
}
```

## References

- Hou, J., et al. (2025). "Efficient estimation in federated causal inference." (In preparation)
- Smucler, E., et al. (2019). "Unifying general and doubly robust cross-fitting."
- Tan, Z. (2020). "Model-assisted inference for treatment effects."

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## Contact

- **Author**: Sinian Zhang
- **Email**: zhan9381@umn.edu
