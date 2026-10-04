# Run on your own data

[Install RoCE](../README.md#install), then copy and edit
[`scripts/real_data/custom.R`](../scripts/real_data/custom.R). The runner reads a
local RDS `data.frame`, including an analysis extract prepared inside an All of
Us workspace. Cohort definitions and confounder selection are dataset-specific.

## Data and preprocessing

Use one row per independent patient and set these fields in `build_data()`:

| Field | Meaning |
|---|---|
| `site_column` | Target/source membership |
| `target_site` | Population whose effect is being estimated |
| `treatment_column` | Numeric 0/1 treatment; contrast is 1 minus 0 |
| `outcome_column` | Numeric outcome: binary (`binomial`) or continuous (`gaussian`) |
| `feature_columns` | Ordered numeric confounder/basis columns shared by all sites |

`site_data_from_frame()` maps the target to `t` and sources to `s1`, `s2`, ...,
preserving row order and recording labels. Outcome and weighting models share
the same features. Exclude IDs, site, treatment, outcome and an intercept column
from the feature list.

Features must be finite and use a common encoding. The template does not impute
or standardize data. Learn any data-dependent transformations on the allowed
training rows and supply their evaluation transforms through `precomputed_folds`.
The [RHC profile](RHC.md) demonstrates outer-fold preprocessing; its inner CV
shares the outer-training transform. The runner handles binary treatment and
uncensored outcomes; it is not a survival-analysis interface.

## Run and configure

```bash
export ROCE_DATA_RDS=/path/to/analysis_data.rds
export ROCE_CORES=2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
Rscript scripts/real_data/run.R /path/to/my_analysis.R /path/to/output/my_run preflight
Rscript scripts/real_data/run.R /path/to/my_analysis.R /path/to/output/my_run all
```

The template selects two-layer, one-round, joint TATE, calibrated nuisances,
`min` CV and 100 lambdas. `lambda_selection=.5` gives aggregation cutoff 2;
it is separate from nuisance penalties. Source/target radii 3/5 reproduce the
RHC profile and need justification for other data.

For three layers, also set `source_validation_method="calibrated"`. Optimizers
are `proximal_newton` or `coordinate_descent`; aggregation also supports
`separate_arms` and `common_tate`. `ROCE_CORES` controls source workers; arms run
sequentially. The implementation emulates federation within one environment.

Replace `all` with `roce`, `sample_size`, `inverse_variance`, `federated_dr`, or
`pooled_dr` to run one method. Use separate output directories for concurrent
processes. Reruns skip completed methods and resume RoCE nuisance checkpoints
with the same input, configuration and build; interrupted baselines restart.
Keep builds, temporary files and results on scratch when using MSI.

## Results

`methods.csv` contains estimates, SEs and 95% CIs on the outcome scale. With the
supplied profiles, `Target-only` is the calibrated target anchor. RoCE also
exports overall/fold-specific weights and nuisance diagnostics; the output
root records configuration, input/build fingerprints, site sizes and R session
information.

Full fits can contain patient-level scores and are saved only with
`ROCE_SAVE_FIT=1`. Outputs stay in your environment; share permitted summaries.
The separate research prototypes in `diagnosis/next_paper/` are not used here.
