# RoCE: Robust Federated Causal Estimation

[![R Version](https://img.shields.io/badge/R-%3E%3D4.0.0-blue.svg)](https://www.r-project.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

RoCE estimates the target average treatment effect from target and source data:

$$\tau_t = E_t\{Y(1)-Y(0)\}.$$

The current analysis uses calibrated nuisance models, two-layer cross-fitting,
and joint optimization of two arm-specific source-weight vectors. Three-layer
cross-fitting, alternative aggregation rules, and one-/two-round communication
are configurable. The R implementation runs site data together as a federated
emulation.

## Install

Requires R >= 4.0 and a C++ compiler.

```bash
git clone https://github.com/h7nian/RoCE.git
cd RoCE
Rscript -e 'install.packages(c("Rcpp", "RcppEigen", "glmnet", "digest", "doParallel"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
```

## Run

Check the installation with synthetic data and all four baselines:

```bash
Rscript scripts/real_data/run.R scripts/real_data/synthetic.R /path/to/output/smoke all
```

For an analysis, follow the [custom-data guide](docs/REAL_DATA.md) for All of Us
or your own dataset, or the [RHC reproduction guide](docs/RHC.md). These use
explicit profiles rather than historical package defaults.

The main R API is `run_tate_crossfit()`. See `?run_tate_crossfit` and the
[configuration template](scripts/real_data/custom.R) for the complete arguments.

| Control | Options |
|---|---|
| Cross-fitting | `crossfit_layers = 2L` or `3L` |
| Aggregation | `joint_tate`, `separate_arms`, `common_tate` |
| Communication | `one_round`, `two_round` |
| Optimizer | `proximal_newton`, `coordinate_descent` |

Cross-fitting layers and communication rounds are separate choices. Joint TATE
retains within-site arm covariance. Baselines include target-only, sample-size
weighting, inverse-variance weighting, Federated-DR, and Pooled-DR.

## Simulations and development

The default simulation DGP uses bounded covariates, four active slopes, and a
fixed population TATE. See the [DGP protocol](diagnosis/repair/bounded_dgp.md)
and [checkpoint guide](diagnosis/repair/checkpoint_recovery.md) for MSI runs.

After installing `testthat` and `withr`, test the installed package from this
checkout:

```bash
ROCE_TEST_INSTALLED=1 Rscript -e 'testthat::test_dir("tests/testthat")'
```

| Directory | Contents |
|---|---|
| `R/`, `src/` | Estimators and C++ solvers |
| `scripts/real_data/` | Portable runner and analysis profiles |
| `tests/` | Regression tests |
| `diagnosis/repair/` | Implementation and validation records |
| `diagnosis/next_paper/` | Separate weak-bias and growing-site research prototypes |

[MIT License](LICENSE) · Sinian Zhang · [zhan9381@umn.edu](mailto:zhan9381@umn.edu)
