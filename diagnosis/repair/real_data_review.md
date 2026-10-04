# RHC application review and current-method alignment

The user requested this review while the MC200 simulations continue. No new
treatment-effect fit was run in this preflight and no historical result was
changed. New artifacts are under scratch `real_data/rhc/`, separate from
simulations and next-paper research.

## Verified input and historical cohort

The vendored RHC CSV has5735 rows and63 columns. Its SHA256 matches the earlier
audit: `9ef4ab578be4b40ad5d97d3a7e08ffdc1f9f76aeeefee51b4996e4221556f8e8`.
The [data provider's description](https://hbiostat.org/data/repo/rhc) identifies
`swang1` as RHC on the first SUPPORT-eligible day and documents the variables.
The [variable catalogue](https://hbiostat.org/data/repo/crhc) lists day1 clinical
measurements; their temporal relation to treatment still needs careful review.

The historical support-screened insurance analysis uses30-day mortality and
Private as target, excluding No insurance and Medicare & Medicaid. This review
reproduced its5039 rows and61 working features:

| Role | Insurance group | n | Treated | Control |
|---|---|---:|---:|---:|
| Target | Private |1698|731|967|
| Source1 | Medicare |1458|511|947|
| Source2 | Private & Medicare |1236|490|746|
| Source3 | Medicaid |647|193|454|

These insurance strata emulate sites; they are not separate hospitals. The
loader's legacy `K=4` means four total groups and three sources. Use explicit
`n_sites` and `source_count` labels in new drivers/reports; do not relabel this
as simulation K4 or force these sample sizes to1000.

Outcome and weight working matrices agree exactly. The pooled intercept-plus-
feature design has rank62; the Medicaid-only design has rank61. There are
2651 missing urine-output values in the retained cohort, currently filled by
the cohort-level median. Seven source/arm/feature combinations are constant
within a source arm but have a different target feature mean. These are finite-
sample moment-support flags, not proof of population nonoverlap. They require
inspection before treating the historical exclusions as sufficient support.

Exact counts and feature flags are in `preflight_20260923_v1/` on scratch.

## Interface differences that require a deliberate update

`run_rhc_tate_experiment()` still calls `run_tate_crossfit()` without explicit
layer, aggregation-mode, target-anchor or score-calibration controls. Its current
defaults therefore invoke the legacy common-TATE/ordinary-target/initial-score
path rather than the new simulation profile. The historical launcher documents
cutoff1, radius5 and its min/conditional-full-1SE numerical protocol.

The result schema also assumes a single source-weight vector: `pairwise$weight`
and `target_anchor_weight` are scalars per source/analysis. Joint TATE needs
separate mu1/mu0 source and target weights. Merely forwarding a joint-mode flag
would produce incorrect labels or recycled rows. A Hou target anchor must also
be distinguished from the ordinary target-only benchmark in method tables.

The current-method application should explicitly record two layers, joint TATE,
score-derivative source calibration, Hou target calibration and cutoff2, with
the existing alternatives retained as named sensitivities. Fitting/inference
radii and nuisance min/1SE rules require the RHC numerical/support audit; the
simulation's radius12 should not be silently substituted for historical radius5.
Preserve the old wrapper defaults and output interpretation for reproducing
archived runs, and test new argument forwarding and arm-specific output schemas.

## Next analyses

1. Review the cohort exclusions, variable timing, missingness handling and
   pooled preprocessing; document what uses covariates versus outcomes and how
   it relates to cross-fitting and federated computation.
2. Inspect the named empirical support gaps, initial density-fit feasibility,
   fitted propensity/transport extremes, effective sample sizes, clipping and
   calibration convergence. Fit audits are separate from favorable effect sizes.
3. Run a checked current-method pilot using the observed unequal site sizes and
   shared feature map. Compare ordinary/calibrated target-only and the four
   existing borrowing baselines on matching records and evaluation folds.
4. Evaluate calibration and aggregation sensitivities, prespecified cutoffs,
   numerical radii, source leave-one-out analyses and fold-seed stability.
   Repeated splits do not create independent cohorts or known-truth coverage.
5. Audit both-arm covariance and eta sensitivity, and distinguish comparison-
   method multiplier draws from full-pipeline refitting. Check reproducibility
   before rendering a new forest plot and source-weight tables.

Real-data comparisons concern point estimates, intervals, borrowing patterns
and robustness. A smaller SE alone does not demonstrate lower bias or valid
coverage. Any causal interpretation remains conditional on observational
identification and transport assumptions. Existing manuscript numbers remain
historical until the corresponding current pipeline is run and checked.

## Current submission, 2026-09-25

The user explicitly authorized real-data jobs. The current wrapper and worker now support Two-layer/joint-TATE, arm-specific weights, distinct ordinary/calibrated target references, and explicit fold seeds. Twelve frozen jobs were submitted under `real_data/rhc/current_two_layer_v1/`; they are not completed results. See [the frozen analysis design](real_data_current_design.md). Old wrapper defaults and old results remain. The final interface gate is `current_interface_checks_20260925_v3/`, with29 interface assertions, existing RHC regressions and an actual synthetic fit/variance decomposition check. The old launcher is retained for historical reproduction; the new worker reuses the checked library and original job-ownership/checkpoint guard.
