# Run RoCE on collaborator data

This entry point accepts a local analysis dataset, including an extract prepared
inside an All of Us workspace or a collaborator's own environment. RHC is one
reproduction example. The current release implements the current-paper RoCE
estimator; `diagnosis/next_paper/` contains separate weak-bias/growing-K research
prototypes and is not used by this runner.

## Install and check

R >= 4.0 and a C++ compiler are required (Rtools on Windows, Xcode command-line
tools on macOS, or the R development toolchain on Linux).

```bash
git clone https://github.com/h7nian/FACE-HD.git
cd FACE-HD
Rscript -e 'install.packages(c("Rcpp", "RcppEigen", "glmnet", "digest", "doParallel"), repos="https://cloud.r-project.org")'
R CMD INSTALL .
# Use a writable output location in your own environment.
Rscript scripts/real_data/run.R scripts/real_data/synthetic.R /path/to/output/smoke all
```

An optional `R CMD INSTALL --library=/path/to/Rlib .` and matching `R_LIBS_USER`
allow installation outside the default library. Generated native registration
files are committed; cloning and installing does not require `devtools` or
`Rcpp::compileAttributes()`. To test the package, additionally install
`testthat` and `withr`. `RCAL` is optional and is not used by this profile.
On MSI keep builds, R libraries, temporary files and outputs on scratch.

## Define the scientific analysis and map the data

Copy `scripts/real_data/custom.R` and edit it for the intended target population,
treatment, outcome and adjustment set. One row represents one independent
patient. The supplied runner supports binary treatment and either a binary
outcome (`family="binomial"`) or a continuous outcome (`"gaussian"`). It does
not convert a censored event-time outcome into a survival analysis.

Save your local analysis `data.frame` as RDS. In `build_data()`, specify:

| Configuration | Meaning |
|---|---|
| `site_column` | Target/source membership; labels may be character or factor |
| `target_site` | The population whose treatment effect is being estimated |
| `treatment_column` | Numeric 0/1 treatment; the contrast is 1 minus 0 |
| `outcome_column` | Numeric outcome; binary outcomes must be 0/1 |
| `feature_columns` | Ordered numeric confounder/basis columns shared by all sites |

`site_data_from_frame()` makes target `t`, then sources `s1`, `s2`, ... in
first-occurrence order and records the label mapping. It preserves row order
and uses the **same feature map** for outcome and merged-weight models. Do not
include patient IDs, site labels, treatment, outcome or an intercept column as
features. No equal-site-size assumption or RHC insurance grouping is imposed.
Each site's input contains `W_outcome`, `Z_site`, `A`, `Y`, and `n` (plus the
matching `X`/`X_dagger` aliases used by some diagnostics).

The simple template expects finite, pre-specified numeric features. It performs
no automatic missing-value imputation, centering, scaling or category selection.
Use a fixed feature dictionary and externally specified transformations. If
transformations must be learned from this dataset, implement them on the
allowed training rows, supplying the corresponding evaluation transforms and
`precomputed_folds` through `build_data()`. The RHC profile illustrates the
outer-fold mechanism; its initial/calibration CV shares each outer-training
transform and is not fully nested preprocessing. Do not fit preprocessing on
the full cohort and describe it as excluding held-out observations.

```bash
export ROCE_DATA_RDS=/path/inside/your/workspace/analysis_data.rds
export ROCE_CORES=2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript scripts/real_data/run.R /path/to/my_analysis.R /path/to/output/my_run preflight
Rscript scripts/real_data/run.R /path/to/my_analysis.R /path/to/output/my_run all
```

All of Us cohort construction, phenotype definitions and confounder selection
remain dataset-specific tasks. This code reads an already prepared local
analysis dataset; it does not download All of Us data. The current R
implementation emulates the federated algorithm using a list of site data in
one process/environment; it is not a distributed service connecting hospitals.

## Fitting choices

The template explicitly selects Two-layer, one-round communication, joint
TATE aggregation, calibrated target/source models, `min` CV and 100 lambdas.
Joint TATE uses **two arm-specific weight vectors** and retains arm covariance.
`lambda_selection=.5` means aggregation cutoff 2; this parameter is distinct
from nuisance penalties selected by CV. Source radii 3 and target radius 5 are
the RHC profile, not a universal recommendation for other datasets.

Use `crossfit_layers=3L` together with `source_validation_method="calibrated"`
for independent inner validation. Use `nuisance_solver="coordinate_descent"`
for the reference solver or `"proximal_newton"` for the accelerated solver.
`aggregation_mode` also accepts `"separate_arms"` and `"common_tate"`.
Specify scientific/tuning choices before comparing effects and retain the
configuration with every result. The runner never uses a known effect to tune
an analysis.

`ROCE_CORES` controls source workers; the portable runner fits arms sequentially.
A method can be assigned to a separate collaborator/job by replacing `all` with
`roce`, `sample_size`, `inverse_variance`, `federated_dr`, or `pooled_dr` and using
a **different output directory for each simultaneous process**. Completed
methods are skipped on rerun. Core RoCE nuisance checkpoints permit resuming
interrupted fits with the same input/configuration/build. Baselines restart if
interrupted. A mismatch requires a new output directory; old results stay intact.

## Outputs and interpretation

Each method directory has `methods.csv` (estimate, SE, 95% CI), warnings and a
completion marker. RoCE additionally records overall and fold-specific weights
and nuisance convergence diagnostics. `Target-only` with the supplied calibrated profiles means the
**calibrated target anchor**. Changing `target_nuisance_method` changes that anchor. Baselines retain their own fitting and variance
protocols; see the RHC guide for that example.

The output root records input/configuration/package fingerprints, site sizes,
feature count, metadata, copied configuration and `sessionInfo()`. Estimates
are on the outcome scale (multiply binary-outcome differences by 100 for
percentage points). There is no known causal truth, bias or coverage measure
for an observational real dataset.

Patient-level fit objects are **not saved by default**. `ROCE_SAVE_FIT=1`
retains them locally when detailed diagnosis is needed. Use the approved
aggregate outputs from the collaborator's data environment when sharing
results; the runner does not upload data or results anywhere.
