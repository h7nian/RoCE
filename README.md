# Robust Federated Causal Estimation (RoCE) through Calibrated High-Dimensional Models

[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

> **Naming convention**: The paper refers to the method as **RoCE**. The R package
> is named `RoCE` (R packages cannot contain hyphens).

## Overview

RoCE implements federated causal estimators for multi-site data. Its primary
estimator fits arm-specific high-dimensional nuisance functions and jointly
optimizes two arm-specific source-weight vectors for the target average
treatment effect (TATE). The current manuscript analysis selects **two-layer
cross-fitting**: outer evaluation and fold-summed calibration. Three layers add
independent inner validation for aggregation weights and remain available as
an explicit comparison. Communication rounds are controlled separately.

**Real-data collaborators:** start with [the general data guide](docs/REAL_DATA.md)
for All of Us or your own data; [the RHC profile](docs/RHC.md) reproduces the
published Private-target example. Both use one portable runner and explicit
preprocessing, tuning, and output settings.

### Default simulation

The default DGP is `bounded_joint_v3`: bounded, population-standardized
truncated-normal features, four active slopes, and a fixed population TATE.
`p` counts working features directly; no quadratic expansion is added.
Source covariate shift and treatment allocation have separate controls.
See the [design and validation protocol](diagnosis/repair/bounded_dgp.md).

The default launcher uses 1000 observations per site, three-level cross-fitting,
and one independent Slurm job per repeat. Builds, results and logs stay under
`/scratch.global/zhan9381/FACE-HD/`. After a checked installation is available:

```bash
bash main.sh --n-sims 200 --p '10 20 50' --max-in-flight 512 --checkpoint
```

The controller first checks one repeat in each cell, then allows up to 512
pending/running repeat jobs (2 CPUs each). `--max-in-flight` accepts 1–2000;
the scheduler determines actual running concurrency. The bounded launcher now
uses `preempt,msismall,agsmall,amdsmall,amd512,ag2tb,msibigmem,saffo-2tb`, with requeue and persistent nuisance checkpoints.
The listed order does not impose a Slurm scheduling priority. Each repeat requests 2 CPUs, 4 GB total memory and 30 minutes.
The lightweight controller requests one hour and checkpoints its submission
ledger before requeueing. For non-preempting runs, choose public partitions
with `--partitions`; `--no-checkpoint` is available there.
See [checkpoint recovery](diagnosis/repair/checkpoint_recovery.md).
Explicit `face` and
`roce` DGP selections remain available for historical comparisons.

### Target Estimand

The primary estimand is the **target average treatment effect**:

$$\tau_t = E_t\{Y(1)-Y(0)\}=\mu^1_t-\mu^0_t.$$

Use `run_tate_crossfit()` for the primary TATE procedure. It forms the
treated and control influence components for joint TATE variance optimization,
retaining the within-site cross-arm covariance. The lower-level
`run_crossfit(A_val = a)` interface remains available for secondary potential-
outcome means $\mu^a_t=E_t\{Y(a)\}$.

### Key Features

- **Calibrated Nuisance Fitting**: Source transport and target calibration with explicit training boundaries
- **Two Communication Protocols**: Two-round and one-round algorithms
- **Two- or Three-Layer Cross-fitting**: Optional independent inner validation of complete calibrated source and target procedures
- **TATE Aggregation**: Common TATE weights, independently optimized arm weights, or jointly optimized arm weights; all retain the applicable within-site arm covariance
- **C++ Acceleration**: Selectable proximal Newton and coordinate descent, with exact fitting-input caches and compact calibration matrices
- **Comparison Methods**: Target-only, sample-size weighted, inverse-variance weighted, Federated-DR, and Pooled-DR estimators with paired-arm variance

### API Stability and Compatibility

- **Primary API**: `run_tate_crossfit()` with `communication_mode = "one_round"` (manuscript default) or `"two_round"`.
- **Arm-specific API**: `run_crossfit(A_val = 0/1)` for secondary potential-outcome means.
- **Reaggregation**: `reaggregate_tate_crossfit()` and its sensitivity grid preserve the fitted source-weight rule, including `quadratic_bias`; unchanged cutoff and inference radius reproduce the original fit.
- **Internal C++ bindings**: functions ending with `_cpp` are internal implementation/testing interfaces and are **not** part of the package's stable public API contract.

Select the current two-layer joint-TATE configuration explicitly:

```r
fit <- run_tate_crossfit(
  data_split, n_folds = 10L, communication_mode = "one_round",
  target_nuisance_method = "hou_calibrated",
  source_validation_method = "outer_fit", crossfit_layers = 2L,
  aggregation_mode = "joint_tate", lambda_selection = 0.5,
  nuisance_solver = "proximal_newton", calibration_layout = "compact"
)
```

`aggregation_mode` can also be `"common_tate"` or `"separate_arms"`, and can be
changed with `reaggregate_tate_crossfit()` using the same nuisance fits.

The development argument `nuisance_cv_certificate = TRUE` enables sufficient
KKT nonconvergence checks in density-ratio CV. It retains the prescribed grid,
CV tolerance and final optimizer. The effective default is `FALSE`; `NULL`
inherits an enclosing call. The control is available in the cross-fitting and
simulation APIs, and as optional argument21 of `main.R`. Cache provenance and
result metadata record it. See the [validation record](diagnosis/repair/density_cv_certificate.md)
before enabling it in a new campaign; existing frozen campaigns are unchanged.
The repeat launcher accepts `--nuisance-cv-certificate` for a library with the
corresponding checked capability marker, including baseline-only workers.

For the two-versus-three-layer comparison, set `crossfit_layers = 2L` or `3L`
in `run_tate_crossfit()`, `run_crossfit()`, `run_single_simulation()`, or
`run_simulation_study()`. Both choices retain an untouched outer evaluation
fold and fold-specific aggregation weights. Two layers reuse the final outer
models to compute weight-learning scores on the outer training sample; three
layers independently validate fully calibrated inner models there. These
training-score definitions differ from the legacy `"initial"` validation
ablation. Layers and communication rounds are separate controls.

The repeat launcher accepts `--crossfit-layers 2` or `3`; it requires a library
that passed the layer checks and records both arm means, their covariance,
TATE variance and fold-specific weights. Keep the DGP, outer folds, calibration,
cutoff, and nuisance grid fixed when comparing layers. Two layers require at
least three original folds, three layers require four, and this comparison
uses the same ten original folds. Fixed-weight and eta-sensitivity variance
outputs must be checked against empirical variance across paired repeats.

The existing `"lasso"`/`"initial"` package defaults remain available for
controlled legacy comparisons. Results record `crossfit_levels` separately
from `communication_mode`. The three-level candidate and its coverage are
still under validation; see [implementation status](diagnosis/repair/README.md).

An experimental score-matched calibration is separately selected with
`calibration_control = list(recipe = "score_derivative",
target_propensity_initialization = "calibrated", target_radius = 12)`.
This radius belongs to the selected bounded-DGP validation profile, which also
uses `M_tau = M_tau_inference = 12`; it is not a universal setting. This option changes the calibration method and still
requires finite-sample and inference validation.

The direct calibration ablation adds
`source_nuisance_method = "standard"` to `calibration_control` and uses
`source_validation_method = "complete"`. It fits an ordinary source OR and
initial merged-weight model on each allowed training subset, omitting final
source score calibration. `"calibrated"` remains the default source program
and a compatible name for complete validation. The standard ablation supports
one-round communication and either target anchor; see the
[factorial comparison design](diagnosis/repair/calibration_ablation_design.md).

For source workers that are recreated between outer folds, optional
`nuisance_cache_dir` reuses fits through a private per-run directory on scratch.
The supplied parent directory must exist and be writable; `use_lambda_cache`
must remain `TRUE`. The default `NULL` uses memory caching. Exact fitting
inputs and model integrity are checked, and the run cleans up its cache.
Disk caching adds overhead on small problems; the full p100 comparison is
still running. Grid reduction is evaluated separately because it can change
the selected penalty and estimator.

## Installation

### Prerequisites

Ensure you have the following R packages installed:

```r
install.packages(c("Rcpp", "RcppEigen", "glmnet", "digest", "doParallel"))

# For comparison methods (optional)
install.packages("RCAL")

# For testing (optional)
install.packages(c("testthat", "withr"))
```

### From Source

```bash
git clone https://github.com/h7nian/RoCE.git FACE-HD
cd FACE-HD
R CMD INSTALL .
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
# partition, retaining their covariance and learning joint arm-specific weights.
result_tate <- run_tate_crossfit(
  data_split,
  communication_mode = "one_round",
  n_folds = 5,
  aggregation_mode = "joint_tate", crossfit_layers = 2L,
  target_nuisance_method = "hou_calibrated",
  source_validation_method = "outer_fit",
  lambda_selection = 0.5, # cutoff c = 2
  verbose = TRUE,
  M_tau = 5,
  M_tau_inference = 5,
  n_cores = 3               # source-site workers per treatment arm
)

cat("TATE estimate:", result_tate$estimate, "\n")
cat("SE:", result_tate$se, "\n")
cat("95% CI:", result_tate$ci_lower, result_tate$ci_upper, "\n")
print(result_tate$weights_by_arm)
cat("True TATE:", data$mu1_true - data$mu0_true, "\n")
```

The current manuscript analysis uses aggregation multiplier `0.5` (Wald
activation cutoff `c = 2`). Historical package defaults remain available;
use explicit arguments or the RHC runner to reproduce the selected analysis. For a pre-specified sensitivity analysis, callers
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

The two-round variant uses each source's initial outcome model. The target
then computes derivative moments for those source-specific models before
source calibration. This introduces the additional communication dependency.

### Optional Three-Layer Cross-fitting

For each outer evaluation fold `k1`, each inner validation fold `k2` evaluates
models calibrated using the remaining folds `k3`. Each calibration block's
initial models exclude `k1`, `k2`, and that block's `k3`. All `k3` losses are
combined into one final model per site, arm, and nuisance type. All `k2`
validation moments are pooled before learning the weights for `k1`.

The final outer model uses calibration on all folds except `k1`, with its own
fold-specific initial models. The three fold roles can use the same ten-fold
partition and do not add communication round trips. See the
[algorithm description](diagnosis/repair/algorithm.md) for sample sizes and
the current inference limitations.

## Grouped nuisance CV for resampling

For resampling diagnostics, a site in `data_split` may optionally carry a
positive integer `cv_group_id` vector identifying original observations.
Repeated origins are assigned together in target and source nuisance CV.
Missing (`NULL`) or all-unique IDs retain the ordinary CV path and RNG behavior.
The supplied outer folds must also keep each origin intact, and their data
references must contain the same group metadata; stale views or cross-fold
origins fail explicitly. For copied data, supply group-consistent
`precomputed_folds` rather than repartitioning copied rows.

This option changes nuisance-validation bookkeeping, not the selected TATE
aggregation objective or its variance formula. It is not a cluster-robust
variance estimator and does not by itself validate bootstrap confidence
intervals. Grouped source CV retains the existing equal-fold score criterion;
groups are balanced by origin count, so row counts per fold can differ.

## Historical FACE Simulation Configurations

The following describes explicit `dgp_type="face"` historical runs. The current
bounded default is described [above](#default-simulation) and in the linked
design protocol.

Both nuisance models use the same working basis $\phi(X)=[X-\kappa, X^2]$ in
every configuration. Misspecification is placed in the true mechanism: the
four signal coordinates are replaced by standardized Kang--Schafer-style
transforms $X^\dagger$ and the true predictor is mixed as
$(1-\omega)\,\eta(X)+\omega\,\eta(X^\dagger)$ with
$\omega=$ `FACE_MISSPECIFICATION_STRENGTH` $=0.75$ (pre-registered, see
`diagnosis/HISTORY.md` #0002):

| Config | True outcome mechanism | True treatment mechanism | Description |
|--------|------------------------|--------------------------|-------------|
| C1 | $\eta(X)$ | $\eta(X)$ | Neither mechanism uses the transformed predictor |
| C2 | $\eta_\omega(X)$ | $\eta(X)$ | Outcome model misspecified |
| C3 | $\eta(X)$ | $\eta_\omega(X)$ | Site/treatment model misspecified |
| C4 | $\eta_\omega(X)$ | $\eta_\omega(X)$ | Both misspecified |

These labels describe the outcome and treatment mechanisms. The actual
propensity also uses the DGP's truncation, and the skew-normal transport ratio
need not belong to the working log-quadratic family. C1 therefore does not
assert that every fitted nuisance is exactly specified.

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

## Historical HPC Campaign (SLURM)

This section preserves the earlier FACE campaign. Use `main.sh` and the
current bounded-DGP protocol for current campaigns.

The historical rerun was restricted to `p=100`. See
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
