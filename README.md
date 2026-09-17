# Robust Federated Causal Estimation (RoCE) through Calibrated High-Dimensional Models

[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Naming convention**: The paper refers to the method as **RoCE**. The R package
> is named `RoCE` (R packages cannot contain hyphens).

## Overview

RoCE implements federated causal estimators for multi-site data. Its primary
estimator fits arm-specific high-dimensional nuisance functions, contrasts the
treated and control influence components, and then learns one common set of
source weights for the target average treatment effect (TATE). The package
provides communication-efficient protocols with two-level cross-fitting.

### Target Estimand

The primary estimand is the **target average treatment effect**:

$$\tau_t = E_t\{Y(1)-Y(0)\}=\mu^1_t-\mu^0_t.$$

Use `run_tate_crossfit()` for the primary TATE procedure. It forms the
treated-minus-control influence values before variance estimation and source
aggregation, retaining the within-site cross-arm covariance. The lower-level
`run_crossfit(A_val = a)` interface remains available for secondary potential-
outcome means $\mu^a_t=E_t\{Y(a)\}$.

### Key Features

- **Doubly Robust Estimators**: Calibrated loss functions for Neyman Orthogonality
- **Two Communication Protocols**: Two-round and one-round algorithms
- **Two-Level Cross-fitting**: Enhanced robustness through nested sample splitting
- **TATE Aggregation**: One common source-weight vector selected from the TATE variance and Wald discrepancies; the smooth quadratic-bias rule (`screening_rule = "quadratic_bias"`) is computed alongside as a pre-specified sensitivity estimator
- **C++ Acceleration**: High-performance coordinate descent with GLMNET-style optimizations
- **Comparison Methods**: Target-only, sample-size weighted, inverse-variance weighted, Federated-DR, and Pooled-DR estimators with paired-arm variance

### API Stability and Compatibility

- **Primary API**: `run_tate_crossfit()` with `communication_mode = "one_round"` (manuscript default) or `"two_round"`.
- **Arm-specific API**: `run_crossfit(A_val = 0/1)` for secondary potential-outcome means.
- **Reaggregation**: `reaggregate_tate_crossfit()` and its sensitivity grid preserve the fitted source-weight rule, including `quadratic_bias`; unchanged cutoff and inference radius reproduce the original fit.
- **Internal C++ bindings**: functions ending with `_cpp` are internal implementation/testing interfaces and are **not** part of the package's stable public API contract.

## Installation

### Prerequisites

Ensure you have the following R packages installed:

```r
install.packages(c("Rcpp", "RcppEigen", "glmnet", "doParallel"))

# For comparison methods (optional)
install.packages("RCAL")

# For testing (optional)
install.packages("testthat")
```

### From Source

```bash
git clone https://github.com/sinianzhang/RoCE.git
cd RoCE
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

# Run the primary TATE estimator. Both treatment arms use the same fold
# partition and are contrasted before learning one source-weight vector.
result_tate <- run_tate_crossfit(
  data_split,
  communication_mode = "one_round",
  n_folds = 5,
  lambda_selection = RoCE:::AGG_WALD_LAMBDA, # selected cutoff c = 1
  verbose = TRUE,
  M_tau = 5,
  M_tau_inference = 5,
  n_cores = 3               # source-site workers per treatment arm
)

cat("TATE estimate:", result_tate$estimate, "\n")
cat("SE:", result_tate$se, "\n")
cat("95% CI:", result_tate$ci_lower, result_tate$ci_upper, "\n")
cat("Common source weights:", result_tate$weights, "\n")
cat("True TATE:", data$mu1_true - data$mu0_true, "\n")
```

The manuscript analysis uses aggregation multiplier `1` (Wald activation
cutoff `c = 1`), selected by the documented coverage-blind rule on disjoint
pilot seeds. For a pre-specified sensitivity analysis, callers
can instead set `lambda_selection = "cv"` and pass a positive
`aggregation_lambda_grid`; candidate weights and the validation criterion are
then computed entirely within each outer-training sample. The grid must not be
supplied together with a fixed numeric `lambda_selection`.

## Project Structure

```
RoCE/
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
│   ├── simulation_diagnostics.R    # Per-setting coverage/RMSE quality control
│   ├── real_data_rhc.R             # RHC preprocessing and analyses
│   ├── numerical_utils.R           # Pure numerical/math functions
│   ├── validation.R                # Input validation functions
│   └── checkpoint.R                # Checkpoint/restart for SLURM preemption
├── src/                            # C++ source files (Rcpp/RcppEigen)
│   ├── density_ratio.cpp           # Density ratio fitting and CV
│   ├── outcome_model.cpp           # Outcome model fitting, CV, and GLM utilities
│   ├── variance.cpp                # Variance, covariance, and correction terms
│   ├── weight_optimization.cpp     # Weight optimization and aggregation
│   ├── optimization.h              # C++ header declarations
│   ├── cv_utils.h                  # Cross-validation shared utilities
│   ├── numerical_constants.h       # Numerical constants (C++ side)
│   └── utils.h                     # GLM utilities and numerical helpers
├── tests/                          # Test suite (testthat)
│   ├── testthat.R                  # Test runner
│   └── testthat/                   # Individual test files
├── scripts/slurm/                  # Bounded p=100 MSI workflows and diagnostics
├── docs/
│   ├── main.tex                    # Mathematical methodology document
│   └── supplemental.tex            # Proofs, assumptions, and supplemental results
├── DESCRIPTION                     # R package description
├── LICENSE                         # MIT License
└── README.md                       # This file
```

## Algorithms

### One-Round Communication Protocol (default)

The target supplies fold-specific initial outcome fits and calibration summaries;
each source returns calibrated nuisance and influence summaries in one exchange.

### Two-Round Communication Protocol

The optional two-round variant lets each source initialize its own outcome
model before target calibration. It relaxes the one-round alignment condition
at the cost of an additional communication exchange.

### Two-Level Cross-fitting

Both protocols use nested sample splitting to provide:
- Neyman Orthogonality through calibrated loss functions
- Cross-calibrated plug-in order where calibrated α uses γ_init as the weight plug-in
- Outer-fold separation between weight learning and evaluation
- A TATE variance computed from within-site centered, treated-minus-control pseudo-values, which automatically includes cross-arm covariance, plus the delta-method contribution of the learned source weights (`se`; the fixed-weight version is returned as `se_fixed_weights`)

## Grouped nuisance CV for resampling

For resampling diagnostics, a site in `data_split` may optionally carry a
positive integer `cv_group_id` vector identifying original observations.
Repeated origins are assigned together in target and source nuisance CV.
Missing (`NULL`) or all-unique IDs retain the ordinary CV path and RNG behavior.
The supplied outer folds must also keep each origin intact, and their data
references must contain the same group metadata; stale views or cross-fold
origins fail explicitly. For copied data, supply group-consistent
`precomputed_folds` rather than repartitioning copied rows.

This option changes nuisance-validation bookkeeping, not the common TATE
aggregation objective or its variance formula. It is not a cluster-robust
variance estimator and does not by itself validate bootstrap confidence
intervals. Grouped source CV retains the existing equal-fold score criterion;
groups are balanced by origin count, so row counts per fold can differ.

## Simulation Configurations

Both nuisance models use the same working basis $\phi(X)=[X-\kappa, X^2]$ in
every configuration. Misspecification is placed in the true mechanism: the
four signal coordinates are replaced by standardized Kang--Schafer-style
transforms $X^\dagger$ and the true predictor is mixed as
$(1-\omega)\,\eta(X)+\omega\,\eta(X^\dagger)$ with
$\omega=$ `FACE_MISSPECIFICATION_STRENGTH` $=0.75$ (pre-registered, see
`diagnosis/HISTORY.md` #0002):

| Config | True outcome mechanism | True treatment mechanism | Description |
|--------|------------------------|--------------------------|-------------|
| C1 | $\eta(X)$ | $\eta(X)$ | Both models correctly specified |
| C2 | $\eta_\omega(X)$ | $\eta(X)$ | Outcome model misspecified |
| C3 | $\eta(X)$ | $\eta_\omega(X)$ | Site/treatment model misspecified |
| C4 | $\eta_\omega(X)$ | $\eta_\omega(X)$ | Both misspecified |

Two deviation mechanisms make source $s_1$ non-transportable by $\rho$ on
the log-odds scale (`deviation_mechanism` in `generate_face_data()` and
`run_single_simulation()`, `HISTORY.md` #0009):

| Mechanism | What moves at $s_1$ | Experiment | Positive-$\rho$ reuse |
|-----------|---------------------|------------|------------------------|
| `treated_arm` (default) | treatment log-odds shift $\Delta_{s_1}=1+\rho$ | `negative_transfer` (C1--C3) | treated arm of $s_1$ refitted |
| `both_arms` | both potential-outcome arms shifted by $\rho$ | `shared_shift` (C1) | both arms of $s_1$ refitted |

Every replicate reports the frozen production row set
`RoCE:::.tate_production_method_rows()`; the quadratic-bias sensitivity row is
gated by `include_quadratic_bias_rule = TRUE`.

## Running Simulations on HPC (SLURM)

The manuscript rerun is restricted to `p=100`. See
[`scripts/slurm/README.md`](scripts/slurm/README.md) for the isolated package
gates, exact same-seed rho-reuse audit, bounded grouped submissions,
aggregation, and per-setting coverage/RMSE diagnostics. The 27,000 canonical
result rows of the negative-transfer family are represented by 4,500
seed/configuration/K jobs, each covering the six rho values; the shared-shift
family adds 9,000 rows in 1,500 jobs under
`results/direct_tate_mc500_b5000/shared_shift/`. Production defaults to one
job at a time and never submits the full grouped manifest in one call. Every setting uses 500 Monte Carlo
replicates; comparison-method intervals use 5,000 paired multiplier-bootstrap
draws. Checkpoints diagnose coverage, RMSE, bias, empirical-versus-reported
standard errors, sparse treatment/outcome cells, density-ratio clipping, and
nuisance/weight-optimizer health separately for every setting and method.

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
@software{roce2026,
  author = {Zhang, Sinian},
  title = {RoCE: Federated Adaptive Causal Estimation in High Dimensions},
  year = {2026},
  url = {https://github.com/sinianzhang/RoCE}
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
