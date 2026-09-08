# RoCE TATE common-weight investigation

## Current review entry points

This directory includes a chronological investigation; older jobs and numerical
results below are retained as historical evidence, not current production.

- [Grouped-CV implementation and current jobs](GROUPED_CV_IMPLEMENTATION.md):
  v19 software gates, actual origin-partition audits and the bounded refit pilot.
- [Stepwise function and intermediate review](STEPWISE_FUNCTION_REVIEW.md):
  input checks, nuisance tuning, influence values and independent reconstruction.
- [Exact aggregation/output contract](AGGREGATION_OUTPUT_CONTRACT.md):
  fold-specific weights, unequal fold sizes and adaptive-variance qualifications.
- [Monte Carlo evidence](N100_WEIGHT_UNCERTAINTY_REVIEW.md): versioned RMSE,
  coverage and fixed-nuisance weight diagnostics.
- [Full completion checklist](COMPLETION_CHECKLIST.md) and
  [RHC preflight](REAL_DATA_FINAL_PREFLIGHT.md): remaining final experiments.

The current grouped-CV package passes the full installed R tests and R CMD
check. Numerical identity and first nonidentity audits do not establish
adaptive-weight inference validity. Final production experiments and real-data
inference remain incomplete.

## Question

Why did the first five C(1), `p=100`, `K=2` production seeds show larger
high-shift RMSE and smaller reported SE than the arm-specific diagnostic and
target-only estimator?

This investigation is diagnostic. It must not tune the estimator on production
coverage or select favorable seeds.

## Estimator contract

For source `j`, RoCE first forms a source-assisted TATE

```text
tau_j = mu1_j - mu0_j,
```

then aggregates the target-only and source-assisted TATEs within each outer
fold, using one common source-weight vector learned on its training complement:

```text
tau_roce,k = (1 - sum_j eta_jk) * tau_target,k + sum_j eta_jk * tau_jk.
```

The final observation-level aggregation uses the actual target/source fold
fractions; with unequal site-specific fractions it is not simply an average
of this display. See the exact finite-sample formula in the
[aggregation contract](AGGREGATION_OUTPUT_CONTRACT.md). Returned mean `weights`
are reporting summaries, not coefficients to apply to pooled `source_estimates`.

There is one common `eta_jk` per source-specific TATE and outer fold. The arm-specific
construction, which learns separate weights for `mu1` and `mu0`, is retained as
a diagnostic comparator and is not the primary estimator contract.

## Variance audit

For fixed weights, the target-site influence is

```text
(1 - sum(eta)) * varphi + sum_j eta_j * zeta_j,
```

and source `j` contributes `eta_j * xi_j`. The implemented quadratic variance,
gradient, Hessian, Wald discrepancy variance, site scaling, and within-site
centering agree with this influence representation. Treated-control covariance
is included by differencing the arm influence vectors before recomputing every
variance and covariance component.

Exact tests cover:

- contrast variance `V1 + V0 - 2 Cov(1,0)`;
- expanded aggregation variance versus direct sitewise pseudo-values;
- `eta = 0` recovery of target-only variance;
- returned variance versus the assembled observation-level pseudo-values;
- foldwise Wald denominator reconstruction.

## Current hypotheses

1. **Structural information loss.** The DGP changes the deviated source's
   treated arm while its control arm remains useful. A common TATE weight must
   downweight both pieces, whereas arm-specific weighting can retain the valid
   control information.
2. **Residual borrowing after detection.** Crossing the Wald cutoff activates
   an L1 penalty but does not guarantee an exact zero weight.
3. **Adaptive-weight uncertainty.** The plug-in influence variance treats
   outer-training weights as conditionally fixed. Finite-sample weight
   selection variability may explain reported SE below empirical SD.
4. **Monte Carlo noise.** Five seeds are insufficient for coverage or stable
   empirical-SD conclusions.

No missing treated-control covariance or sign/scaling error has been found in
the inspected fixed-weight algebra. This does not exclude a first-order
adaptive-weight contribution to the unconditional sampling distribution.

## Required intermediate diagnostics

Every new RoCE TATE row records, by source and outer fold:

- common TATE weight;
- Wald statistic and penalty coefficient;
- source-target TATE discrepancy;
- fold-to-fold weight variability;
- corresponding `mu1` and `mu0` arm-specific weights and diagnostics.

The aggregation script is
`scripts/slurm/summarize_tate_source_diagnostics.R`.

## Staged experiment gates

Use one `(config, K, seed)` per Slurm array task. Each task reuses the same
generated data and reference nuisance fits across all six rho values and commits
the full bundle before writing a sentinel.

1. One-seed smoke: verify complete bundle, provenance, finite diagnostics, and
   source identities.
2. Ten seeds: inspect source-specific filtering and paired common/arm-specific/
   target-only errors. This is diagnostic, not confirmatory.
3. Twenty-five and fifty seeds: inspect paired MSE, bias/MCSE, empirical SD,
   mean SE, coverage MC interval, and source-weight behavior.
4. If the SE ratio remains materially below one or common weighting remains
   worse at 50 seeds, pause expansion and run frozen-weight/oracle-source-mask
   experiments.
5. Only after the decision gate passes, expand through C(1) K=4/8 and then
   C(2)--C(3), ultimately reaching 500 seeds per setting.

## Historical initial job DAG (not the current queue)

- `18037057_1`: C(1), K=2, seed 1, six-rho smoke.
- `18038561_[2-10%2]`: seeds 2--10, `afterok:18037057`, at most two concurrent.
- `18038628`: source-diagnostic aggregation, `afterok:18038561`.
- `18038719`: six-rho checkpoint audit, `afterok:18038561`.

Outputs are isolated under
`results/direct_tate_mc500_b5000/roce_tate_diagnostic_v1/` and must not be
combined with the older v35 production fingerprint.

## Seed-1 smoke evidence

Job `18037057_1` completed successfully in 1:04:28 and atomically committed all
six rho outputs. Source-specific diagnostics show the intended filtering:

| rho | common weight s1 | common weight s2 | mean Wald s1 | mean Wald s2 |
|---:|---:|---:|---:|---:|
| 0.0 | 0.324 | 0.221 | 0.96 | 0.74 |
| 0.5 | 0.190 | 0.267 | 2.82 | 0.74 |
| 1.0 | 0.041 | 0.316 | 4.68 | 0.74 |
| 1.5 | 0.000 | 0.329 | 5.99 | 0.74 |
| 2.0 | 0.000 | 0.329 | 7.13 | 0.74 |
| 2.5 | 0.000 | 0.329 | 7.78 | 0.74 |

The deviated source is therefore fully removed by rho 1.5, while the
informative source and its diagnostics remain rho-invariant. Estimate and SE
are bitwise identical for rho 1.5, 2, and 2.5. This rules out failed source
detection as the explanation for the seed-1 high-shift behavior.

The smoke also exposed a diagnostics bug: arm-specific aggregation weights did
not carry source names, so `mu1_source_*` and `mu0_source_*` columns were absent.
The dependent seed array was cancelled after two tasks had run for about five
minutes; neither left partial outputs. Their stale locks were quarantined under
`recovery/cancelled_18038561/`. The source-name contract is now fixed and must
pass a fresh package gate before resubmission.

## Seed-1 v2 arm-specific evidence

After preserving source names on arm-specific aggregation objects, the v2
smoke (`18052947_1`) completed in 1:04:18. Every rho output contains finite
common, `mu1`, and `mu0` source/fold diagnostics. At rho 1.5 and above:

| weight system | s1 weight | s2 weight |
|---|---:|---:|
| common TATE | 0.000 | 0.329 |
| arm-specific mu1 | 0.000 | 0.369 |
| arm-specific mu0 | 0.317 | 0.194 |

This directly confirms the structural-information hypothesis for this seed.
The common estimator correctly rejects the deviated source-specific TATE. The
arm-specific diagnostic can additionally retain that source's unchanged
control-arm information, which is unavailable under the primary one-weight-per-
source-TATE contract. This is an efficiency comparison, not evidence that the
common TATE screening failed.

## Ten-seed soft-penalty checkpoint

The v2 array completed ten independent C(1), `p=100`, `K=2` seeds for all six
rho values. All implementation checks passed; the gate intentionally does not
waive statistical-review flags.

| rho | RoCE RMSE | target RMSE | RoCE bias | coverage | mean SE / empirical SD |
|---:|---:|---:|---:|---:|---:|
| 0.0 | 0.0156 | 0.0328 | -0.0035 | 1.00 | 1.38 |
| 0.5 | 0.0301 | 0.0328 | 0.0240 | 0.90 | 1.14 |
| 1.0 | 0.0417 | 0.0328 | 0.0298 | 0.80 | 0.76 |
| 1.5 | 0.0403 | 0.0328 | 0.0221 | 0.80 | 0.70 |
| 2.0 | 0.0371 | 0.0328 | 0.0154 | 0.80 | 0.72 |
| 2.5 | 0.0357 | 0.0328 | 0.0126 | 0.80 | 0.74 |

The deviated source's mean common weight declines from 0.303 at rho 0 to
0.029 at rho 2.5, but only half the rho-2.5 seeds set it exactly to zero. An
approximate post-hoc mask of that source reduces high-rho RMSE to about 0.026,
below target-only, while fixing the informative-source weight adds almost no
further improvement. Because the v2 CSVs lacked exact source estimates, the
post-hoc reconstruction error can be as large as 0.015 and is not sufficient
to change the primary method.

The next diagnostic therefore records exact source estimates and foldwise
discrepancies and evaluates a predeclared hard-Wald rule: sources above the
foldwise cutoff receive zero weight, and variance minimization is performed
only over retained sources. This rule is opt-in and must beat the soft rule in
paired 25/50-seed gates without unacceptable rho-0 loss before adoption.
