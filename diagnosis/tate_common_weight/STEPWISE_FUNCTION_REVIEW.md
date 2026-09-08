# Stepwise function, parameter and intermediate-value review

This review responds to the request to test individual functions and inspect
intermediate values. It covers the critical TATE path on the saved seed-1,
rho-0 identity fit and first genuine resampling draw. It is **not** a claim
that every function, parameter combination, simulation setting, or confidence
interval in the repository has been validated.

Latest version: v19 full installed tests and R CMD check passed. Both grouped-CV
identity refits and the first actual rho=0 grouped draw passed numerical audits.
Its 145 explicit nuisance-CV records contain no cross-fold duplicate origins;
the independently reconstructed common weights, Wald statistics, point and
site-centered variance pass, with maximum normalized KKT residual 4.24e-8.
The four paired plug-in SEs are 0.02120638, 0.02129993, 0.02063135 and 0.02067566.
See `GROUPED_CV_IMPLEMENTATION.md` for the current result and artifact paths.
The detailed table and instability trace below describe the **archived row-CV
draw**, not the corrected grouped-CV result. Neither is a coverage validation.

### Rho-specific audit reference selection

The first grouped rho=1 draw passes the general intermediate auditor and an
independent actual-origin audit (145 valid explicit records, zero cross-fold
origins). Its point reconstruction error is zero and variance error 1.08e-19.
The separate optimizer auditor exposed its own hardcoded `artifacts[["0"]]`
reference selector: applying that script to rho=1 incorrectly compared the
identity fit with rho=0 and stopped before publishing any canonical audit.
Using the actual rho=1 reference gives identity error zero. This is an audit
reference-selection defect, not evidence of a failed estimator or optimizer.
The auditor was corrected and passed 18 boundary/selection tests and fresh
rho=0/rho=1 saved-fit regressions. It uses actual rho, verifies paired metadata
and CSV/RDS agreement, and rejects fractional IDs and malformed scalars.
The rho=1 batch expanded only after these gates passed; no failure threshold
was weakened and no frozen worker changed. Canonical optimizer audits now use
the fresh `optimizer_v2` directories listed in `GROUPED_CV_IMPLEMENTATION.md`.

### Repeated same-resample check after grouping correction

The grouped rho=0 positive draws 1--4 now pass the general intermediate audits.
Independent variance traces for draws 2--4 verify exactly matching row maps and
multiplicities against their archived row-CV counterparts. Their refitted,
relearned-weight plug-in SEs are 0.02063838, 0.02207229 and 0.02176139, versus
0.06853161, 0.03849082 and 0.03851743 under row CV. The target treated-arm
extreme AIPW values shrink as well; see the paired table and artifact paths in
`GROUPED_CV_IMPLEMENTATION.md`. This is reproducible intermediate evidence of
protocol sensitivity across four resamples, not independent Monte Carlo
coverage evidence. The predeclared 20-draw batches continue unchanged.

### Complete rho=0 pilot: retained high-leverage source observation

All 20 rho=0 draws now pass the numerical auditor, but draw 13 has materially
larger within-draw variance. Saved-model/row-map reconstruction locates its
main source-1 tail at original row 776 (A=Y=1), sampled twice at positions 322
and 658 in outer fold 5. Both arm/fold assignments agree; differing local fold
vector positions are not observation mismatches.

For this treated observation, the source correction is
`exp(-g) * (Y - m_alpha)`. The paired factors are:

| Factor | Fixed nuisance | Refit grouped-CV nuisance |
| --- | ---: | ---: |
| Density logit g | -2.907991 | -4.877261 |
| Density ratio exp(-g) | 18.31996 | 131.2707 |
| Outcome prediction | 0.691219 | 0.708832 |
| Outcome residual | 0.308781 | 0.291168 |
| Raw source correction | 5.656851 | 38.22188 |
| Fold-5 source-1 weight | 0.350079 | 0.369412 |
| Site-scaled raw influence, 3 * weight * correction | 5.941029 | 42.35885 |

The source treatment indicator is one and N/n_source remains three. Raw
corrections and phase-3 pseudo-values reproduce their saved values exactly.
The density-ratio factor increases about 7.17-fold while the residual decreases
and the common weight increases only about 5.5%. Thus density tilting, not the
outcome fit or row/arm indexing, drives this tail. The logit is still inside
M_tau=5, just below the permitted ratio exp(5), so absence of a clipping event
does not imply low leverage. No new clipping or fitting rule was introduced
to suppress this draw. See its versioned variance trace under
`full_refit_grouped_cv_v19_audits/seed_000001_rho_0_draw_0013_variance_trace_v1/`.

### Exact grouped-CV warning replay

Before interpreting the historical inventory below, note that the later
draw-12 grouped-CV warning has been localized by exact target-only replay:
75/75 calls reproduce saved predictions and group/fold partitions exactly.
Only target PS at k1=5, k2=NULL emits glmnet -80; its selected lambda.min is
0.0495189105881761. This does not prove absence of tuning-search truncation
effects. The first exploratory replay used the wrong OR family; only the
corrected Binomial `draw12_target_glmnet_warning_replay_v2/` is accepted as
replay evidence. See `GROUPED_CV_IMPLEMENTATION.md` for scope and hashes.

### Historical row-CV component test inventory

| Layer / functions | Checks actually performed | Evidence / status |
| --- | --- | --- |
| `roce_refit_resample` and observation partitions | Integer seeds, finite exact site schema, all fields share row maps, complete fold partitions, RNG preservation, invalid maps/counts and cross-fold copies rejected | 143 resampling/saved-fit assertions plus the separate auditor negative tests; actual draw has 1,878 unique origins in 3,000 positions, max multiplicity 5, zero cross-fold origins |
| Target-only AIPW components | Both arms' saved propensity and outcome predictions reconstructed on original/resampled row IDs; manually recomputed `m + I(A=a)*(Y-m)/p_a` with the package probability floor | All 4,000 original/refitted arm-level target pseudo-values reconstruct; extreme fitted predictions identified below |
| TATE contrast helpers | Shared treated/control outer and inner IDs; estimate and influence contrasts; raw variance/covariance moments | Outer and inner moment maximum discrepancy about 1.25e-16; common-arm reconstruction passes |
| Training Wald statistic and penalty | Training-site sizes, discrepancy variance/SE, Wald values, lambda=1 truncated penalty | Independent reconstruction agrees with stored values |
| Weight optimizer | Normalized quadratic objective, gradient/subgradient and KKT in signed unconstrained R^K; independent optimizer rerun and ridge reconstruction | Max KKT residual 9.08e-7 (<1e-5); rerun weights identical; ridge zero in the two inspected fits |
| `.compute_phase3_all_phi` | Actual fold-specific common weights, N/n_site scaling, each observation assigned once; both-arm subtraction at the same weights | Final point and pseudo-value reconstruction agree; mean fold weights are not misused as global pooled weights |
| Site-centered variance | Within-site centering; independently computed sums of squares; per-site `V1 + V0 - 2 Cov(1,0)` | Final variance errors zero; largest displayed arm-covariance identity error about 2.2e-19 |
| Saved-result audit interfaces | Valid rehashed wrong maps, counts, scalar disagreement, altered moments, invalid dimensions/source IDs, overwrite protection and wrong identity weights | 8 negative/positive test blocks, 13 explicit expectations passed |

The full installed-package test suite and R CMD check passed for the unchanged
v18 package. Those software gates do not replace this numerical review or
establish frequentist coverage.

## A concrete intermediate instability, not an aggregation mismatch

The first resampling draw produced these four **within-draw** plug-in SEs,
all on the same resampled observations:

| Nuisance fits | Common weights | TATE estimate | Plug-in SE |
| --- | --- | ---: | ---: |
| Original, fixed | Original | 0.20850 | 0.02121 |
| Original, fixed | Relearned | 0.20731 | 0.02130 |
| Refit | Original | 0.16852 | 0.04752 |
| Refit | Relearned | 0.18047 | 0.04299 |

These are not SDs across independent bootstrap draws. The comparison localizes
the large plug-in SE change to nuisance refitting in this particular draw;
relearning weights then reduces that refitted plug-in variance. It does not
prove a positive missing variance term or a systematic bootstrap bias.

For target original row 295 (resampled position 773, outer fold 4, A=1,Y=0):

| Intermediate | Original nuisance | Refitted nuisance |
| --- | ---: | ---: |
| Estimated P(A=1|X) | 0.1187423 | 0.0135263 |
| Predicted treated outcome | 0.6605867 | 0.8768983 |
| Treated AIPW pseudo-outcome | -4.90261 | -63.95203 |

This observation is not at the stored propensity clipping boundary. Its large
residual divided by a small fitted propensity explains a large target influence
value. Under refitted nuisances/original weights, one observation contributes
about 50.7% of the target-site plug-in variance; with relearned weights this is
about 45.1% of the smaller target-site contribution. The relearned estimator's
total variance is distributed approximately 25.5% target, 62.9% source 1, and
11.6% source 2, so the target outlier is not the sole contributor.

The next target-CV audit will reproduce just the affected complement-fold fits
and inspect the actual selected penalties and model complexity. Nuisance-CV
duplicate-origin leakage remains a possible mechanism, not an established sole
cause. No clipping threshold, nuisance rule, or primary weight objective has
been changed to make this observation disappear.

## Target-CV replay and controlled grouped-origin check

The targeted fold-4 replay subsequently reproduced the original and resampled
PS, outcome predictions and full AIPW vectors exactly (maximum AIPW rounding
error 2.22e-16). The resample selects weaker penalties and denser models:

| Model | Original lambda.min | Resampled lambda.min | Original/resampled nonzero slopes |
| --- | ---: | ---: | ---: |
| Target PS | 0.0288281 | 0.0128664 | 25 / 95 |
| Target treated OR | 0.0430549 | 0.0276295 | 18 / 51 |

Using glmnet's actual returned CV fold IDs (after verifying `keep=TRUE` does
not alter tuning or predictions here) establishes **actual**, not merely
possible, duplicate-origin overlap in this resampled training set:

- PS: 800 rows, 496 unique origins; 187 origins cross CV folds, involving
  460 rows with copies in other folds.
- OR: 292 rows, 190 unique origins; 65 origins cross CV folds, involving
  156 rows with copies in other folds.

This is nuisance-internal CV overlap in the bootstrap multiset. It is not an
outer or secondary cross-fitting leakage: those partitions passed their exact
ID tests, and original unresampled rows have no such duplicate-origin problem.

On the **same resampled training rows**, assigning nuisance-CV folds by origin
removes all duplicate-origin overlap. In this controlled sensitivity, PS
lambda.min becomes 0.0625643 with 4 nonzero slopes and OR lambda.min becomes
0.0926032 with 2 slopes. The identified row's propensity becomes 0.2050621,
its treated outcome prediction 0.6970449, and its AIPW value -2.70214.
Grouping also changes CV fold composition and effective training partitions;
this does not prove overlap is the sole cause or establish valid bootstrap
inference. It does demonstrate that the naive nuisance-refit reference is
highly sensitive to its CV partition semantics and needs correction/validation
before it can support the primary variance decision.

These additional checks are in `audit_target_fold4_tuning.R` and
`audit_target_fold4_grouped_cv.R`; their immutable outputs are under
`full_refit_calibration_v18_audits/seed_000001_rho_0_draw_0001_target_fold4`
and `..._target_fold4_grouped_cv_v2` respectively. No primary package or
scientific fitting parameter was changed for these tests.

## Read nuisance diagnostics literally

Both inspected fits report zero final nonconvergence, line-search failures,
outcome degeneracies and support-floor activations. Source nuisance coefficient
maxima increase in the resample, but are far below the recorded hard boundary.
Large CV `invalid_fold_fits` counts include deliberately skipped path tails:
the original/actual-draw totals are 16,925/17,425, of which 16,319/16,825 are
tail skips. They are candidate-path diagnostics, not thousands of failed final
models. The CV aggregator rejects candidates without valid losses in all folds;
it errors if no candidate is usable. Selected-loss fields are not present in
the saved 10-by-32 nuisance summary and must not be invented.

## Reproducible artifacts and remaining tests

- Structural/raw-moment audits:
  `results/direct_tate_mc500_b5000/full_refit_calibration_v18_audits/seed_000001_rho_0_draw_000{0,1}_v2/`.
- Four-stage variance and per-observation target prediction trace:
  `results/direct_tate_mc500_b5000/full_refit_calibration_v18_audits/seed_000001_rho_0_draw_0001_variance_trace_v2/`.
- Independent optimizer/Wald/contrast review:
  `diagnosis/tate_common_weight/audit_results/seed1_rho0_draw0_draw1_intermediates_v3/`.
- Executable checks are `test_full_refit_resampling.R`,
  `test_audit_full_refit_intermediates.R`, `audit_full_refit_intermediates.R`,
  `audit_full_refit_aggregation_intermediates.R`, and
  `trace_refit_variance_components.R` in this directory.

Further actual draws and rho=1 identity/resampling checks remain necessary.
The computationally passed rho=0 draw permits the predeclared bounded draw
batch, not the full production experiment or a claim of validated inference.

## Subsequent implementation checkpoint

Before interpreting the grouped protocol, the cached row-CV draws 2--4 were
also traced with the pinned v18 package, without any additional fitting.
All point/variance/arm-covariance identities and output hashes passed. The
within-draw SE pattern seen in draw 1 recurs:

| Legacy draw | Fixed nuisance/original weights | Fixed nuisance/relearned weights | Refit nuisance/original weights | Refit nuisance/relearned weights |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0.02207 | 0.02234 | 0.06891 | 0.06853 |
| 3 | 0.02143 | 0.02132 | 0.03739 | 0.03849 |
| 4 | 0.02128 | 0.02166 | 0.03814 | 0.03852 |

Draw 2's maximum absolute treated target AIPW value is 76.96, with seven
target positions at the stored lower propensity bound; draws 3 and 4 have
maxima 16.17 and 18.24 and three/two lower-bound positions. The PS values are
repeated in the separate arm tables, so those counts must not be doubled.
All three fits have zero recorded final nonconvergence, line-search failures,
degeneracies and support-floor activations. Numerical convergence therefore
does not rule out these influential prediction tails.

These remain old row-level-CV diagnostics. Their four-draw empirical SDs and
covariance decomposition are too imprecise, and the protocol too problematic,
to validate primary inference. The new trace outputs are the corresponding
`seed_000001_rho_0_draw_000{2,3,4}_variance_trace_v1` directories under
`results/direct_tate_mc500_b5000/full_refit_calibration_v18_audits/`.

Following the confirmed internal-CV overlap, opt-in grouped CV was implemented
across target glmnet and all source C++ paths. Full installed tests and
`R CMD check` passed for the new immutable v19 snapshot; ordinary/all-unique
inputs retain their tested numerical behavior. See
`GROUPED_CV_IMPLEMENTATION.md` for exact package hashes, scope and remaining
production-dimensional resampling checks. The old row-CV pending draws remain
held; passing software tests does not make their statistical interpretation valid.
