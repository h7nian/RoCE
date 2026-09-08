# Grouped nuisance-CV implementation and validation gates

## Purpose and scope

Latest numerical checkpoint: both v19 identity gates passed. The first actual
rho=0 draw (`18212040`) completed 0:0 in 00:32:11 and passed the current
intermediate auditor plus independent origin-partition and optimizer reviews.
The first actual rho=1 draw (`18213685`) has also completed and passed the
current intermediate auditor and independent origin/optimizer reviews. The chronology
below preserves earlier states; neither identity nor one resample establishes
coverage or a production-ready variance policy.

Both predeclared positive-draw sets are complete and audited: rho=0 20/20 and
rho=1 20/20, plus two identity draws. The independent completion audit verified
all 42 exact bundles and canonical audit hashes, with no missing/duplicate IDs
or locks. All 6,090 actual CV records are valid: 5,800 explicit positive-draw
records and 290 all-unique identity records. Point reconstruction error is zero,
maximum variance error 1.08e-19 and raw-moment error 3.47e-17. There are no final
nuisance convergence/line-search/degeneracy/support-floor failures and only the
two already investigated draw-12 warnings. This completes the bounded pilot,
not the full goal or inference validation. The next independent-dataset gate
is predeclared in `INDEPENDENT_INFERENCE_PILOT.md`.

Draw 20 completed 0:0 in 00:42:22; its point estimate is 0.1669040864497 and
variance 0.000656568547579733, with zero point/variance reconstruction error,
145 valid explicit partitions and no new warnings. Earlier completion records
are retained below. Draw 19 completed 0:0 in 00:49:25 and
passed its canonical audit/checksums: 145 valid explicit CV records, zero
point/variance error and no final nuisance failures or warnings. Its estimate
is 0.222876315381653 and variance 0.00047100646256468.
Draw 18 completed 0:0 in
00:47:24 and passed its canonical audit and checksums: 145 valid explicit CV
records, zero point/variance error and no final nuisance failures or warnings.
Its estimate is 0.202119329534726 and variance 0.0006278427860701.
After the two completed
draws released their slots, read-only `scontrol` confirmed draw 20 RUNNING
with throttle 4; it was not resubmitted.
Draws 16 and 17 completed 0:0 in 00:42:25 and 00:40:23; both passed canonical
intermediate audits and payload hashes, each with 145 valid explicit CV
records, zero point/variance reconstruction error and no final nuisance
failures or warnings. Their estimates are 0.201779288476724 and
0.192914547502068, with variances 0.00059496804557657 and 0.000556957625357659.
Rho=0 draw 5 (`18217670_5`) completed 0:0 in 00:40:24; its canonical
`seed_000001_rho_0_draw_0005/` audit has 145 valid explicit CV records, zero
point/variance reconstruction error, maximum raw-moment error 6.94e-18, and
zero final nuisance failures or warnings. Its estimate is 0.18883066501754 and
variance 0.000443454263483405; all eight audit payload checksums pass. These
counts exclude identities and do not count still-running draws as completed.

The first row-resampled nuisance-refit diagnostic exposed actual duplicate
origins in both training and validation of target nuisance CV. The working
package now has opt-in group metadata support across the complete TATE
nuisance pipeline. The installed v18 package and its completed scientific
outputs remain unchanged. This is a new development version, not a passed
inference procedure.

## Interface contract

- A site may supply `cv_group_id`, a positive integer vector aligned with all
  its rows. It represents observation origins, not treatment-arm labels or
  site weights. Matrices, fractions, missing/nonfinite values and invalid
  lengths are rejected.
- `materialize_fold()` and `combine_folds()` preserve those IDs through
  target/source training and calibration stacks. Stale precomputed data views
  or an origin spanning outer folds are rejected before fitting.
- Target propensity CV uses all target training IDs; target outcome CV uses
  the appropriate treatment-arm subset. Initial outcome, initial density
  ratio, calibrated density ratio and calibrated outcome models receive the
  matching metadata. Local lambda reuse otherwise remains the existing policy.
- Group-aware target caches are distinct from ordinary caches. With no group
  metadata, existing cache keys are unchanged.
- The shared R group-fold builder keeps every repeated origin in a single
  fold. The five C++ CV selectors accept optional one-based `cv_fold_id`
  values aligned to their **arm-filtered** rows and validate them strictly.
  Explicit partitions return their actual fold IDs and validation sizes.
- `NULL` and all-unique groups consume no additional RNG and use the original
  nuisance CV routines. This behavior is tested, including complete TATE
  point estimates, SEs, common weights and pseudo-values.

Groups are balanced by the number of origins, not by replicate multiplicity.
No hidden change was made to source CV's existing equal-fold scoring/one-SE
criterion. The grouped validation folds can have unequal row counts; this
design remains a resampling protocol requiring numerical and statistical
assessment, not a new claim of inferential validity.

## Checks and current status

Completed in the development checkout:

- R parsing, Rcpp interface regeneration, development compilation and Rd
  regeneration.
- Strict group input tests, NULL/unique RNG neutrality, intact-group tests,
  target PS/OR wiring and cache-isolation tests.
- Full one-round and two-round propagation tests that observe the group-fold
  builder at every target/source tuning stage.
- Complete TATE numerical equivalence for ordinary versus all-unique metadata
  on the same fixture.
- Negative outer-view tests: missing/mismatched metadata, invalid partitions,
  cross-fold origins and invalid group dimensions.

The initial propagation test had an incorrect expected diagnostic label: the
actual target calls distinguish `estimate_complement_fold_aipw PS` and `OR`.
After correcting that expectation, the targeted R suites passed. No numerical
implementation was changed to force that test to pass.

## Immutable v19 software gates

The independent R integration review found no missing forwarding, cache mixing,
or default-path RNG changes. All five C++ selector tests passed (65 assertions).
The final combined targeted test run passed after adding explicit rejection of
complex-valued group IDs as well.

- Installed full-test job `18210135`: COMPLETED 0:0 in 00:03:55; all tests
  passed, with recorded sparse-category warnings retained in its log.
- R CMD check job `18210138`: COMPLETED 0:0 in 00:04:16, exact `Status: OK`.
- Source fingerprint:
  `8e17115d87a5cfa0ba9f4051c4018989eb9c3b96067d5f9eeece792edd4b82ac`.
- Installed package fingerprint:
  `2a6ba02daaadc448e63563bd78eab574a1d80c8edfcfdfa0940a7dcb76f91ec7`.
- Test-suite fingerprint:
  `b262df4a6e6191d47adc8205f4d21689f72b84f071827170d541f7a99b2dc494`.

Installation and gates are under
`results/direct_tate_mc500_b5000/Rlib_grouped_cv_20260905_v19/` and
`package_check_grouped_cv_20260905_v19/`. The source/package hashes were
recomputed after completion and match both gates. These certify the software
snapshot, not a corrected production-dimensional bootstrap distribution.

Remaining numerical/statistical gates (the corrected runner is now implemented
and tested as described below):

1. Identity checks on the original production-dimensional datasets and
   nonidentity audits showing that both target and source nuisance CV honor
   the groups, followed by the predeclared bounded resampling checkpoint.
2. Adequate independent-sampling coverage and RMSE evaluation before changing
   primary inference or launching the full production experiment.

## Corrected diagnostic workflow checkpoint

The current `full_refit_resampling.R` now has an explicit `origin` nuisance-CV
mode. It attaches the original observation IDs to every resampled site and
updates each precomputed fold's data reference before fitting. It records every
call to the installed group-fold builder in an isolated sequential R process,
restores the function on success/error, and stores the actual group/fold vectors
in `draw.rds`. A missing or invalid partition record cannot pass the draw gate.

Standalone validation passed 32 expectations: malformed records, RNG/return
preservation, trace cleanup, preexisting-trace refusal, and small Gaussian
identity/nonidentity refits covering all six nuisance stage labels. The
dual-protocol auditor passed 20 expectations; grouped/mixed-protocol summary
fixtures also passed. The readers recompute origin/fold consistency, not merely
trust the saved `valid` flag. Ordinary archived row-CV results remain explicitly
marked `row_level` and cannot be mixed with `origin_grouped` results.

The controller pins the tested v19 refit package while separately recording the
v18 source-fit fingerprint. Its current ordered workflow fingerprint is
`555fbf44d6966b843305d3c92cf7ce65a50e172624836efd49ae02354036217d`.
This is distinct from the archived row-CV protocol.

Two production-dimensional identity checks were submitted and verified RUNNING:

- `18210893`: seed 1, rho 0, draw 0.
- `18210910`: seed 1, rho 1, draw 0.

Both use five CPUs, the v19 installation, and new destinations under
`results/direct_tate_mc500_b5000/full_refit_grouped_cv_v19/`.

The rho=0 job `18210893` subsequently completed 0:0 in 01:28:28 and passed
both the current numerical auditor and an independent read-only audit. All
four estimates, fitted SE/variance, fold weights and Wald values match the v18
source fit exactly. The 145 actual CV records cover initial outcome (10),
initial density (20), calibrated density (20), calibrated outcome (20), target
PS (25) and target OR (50). Every record is valid; all fold-ID vectors are NULL
because every original-ID vector is unique in this identity draw. Thus the
legacy numerical path is correctly preserved without masking duplicated IDs.
The draw/source hashes and all eight audit payload hashes passed. Its audit
is under `full_refit_grouped_cv_v19_audits/seed_000001_rho_0_draw_0000/`.

After that prerequisite passed, the first actual rho=0 origin-grouped draw
was submitted as job `18212040` (seed 1, draw 1) and verified RUNNING. The
rho=1 identity job `18210910` remains running. The next gate is numerical and
explicit-group-partition auditing of the nonidentity output; the identity
pass is not a bootstrap variance or coverage conclusion.

The rho=1 identity job subsequently completed 0:0 in 01:50:34 and passed the
same current auditor plus independent review. Its fitted estimate, SE,
variance, fold weights and Wald values match the v18 baseline bit-for-bit.
The two matched fixed-nuisance estimates differ by one ULP (2.78e-17), well
below the 1e-10 gate; they must not be described as bit-exact. All 145 CV
records are valid with unique origin IDs and NULL fold overrides, and all
draw/source/audit hashes pass. The canonical audit is
`full_refit_grouped_cv_v19_audits/seed_000001_rho_1_draw_0000/`.

After this gate, the first real rho=1 grouped draw was submitted as job
`18213685` (seed 1, draw 1). It uses the same deterministic resampling seed as
rho=0 draw 1 (`18212040`) to retain paired row maps. Both positive-draw
results still require actual origin-fold and numerical audits before the
predeclared batch is expanded. These identity gates do not establish coverage.

The old row-CV array's draws 2--4 completed and were retained as diagnostic
evidence. Its held draws 5--20 were subsequently cancelled only after verifying
all sixteen remained unstarted; Slurm elapsed times are all 00:00:00. No output is
deleted or overwritten, and no old result counts as corrected grouped-CV
inference evidence.

Before changing the diagnostic runner/helper, the exact prior row-CV protocol
was preserved read-only at
`results/direct_tate_mc500_b5000/full_refit_row_cv_v18_protocol_source/`.
Its 13 payload hashes passed and its ordered workflow hash matches the
recorded v18 row-CV hash `415e970b050460e5b9e35f721c172e0bdef393e4101ac997e70467556e6aa5ed`.
Do not overwrite that archive, the installed v18 library, or earlier outputs.

## First corrected nonidentity draw: seed 1, rho 0, draw 1

All 145 nuisance-CV calls have explicit fold vectors with intact origin groups;
the independent validator found zero origins crossing nuisance-CV folds. Row
maps and multiplicities are identical to the archived row-CV draw, permitting
a paired protocol comparison. Final nuisance fits have zero nonconvergence,
line-search failure, degeneracy, support-floor events, or recorded warnings.
Skipped CV candidate tails remain reported separately from final-fit failures.

| Nuisances | Common weights | Estimate | Within-draw plug-in SE |
| --- | --- | ---: | ---: |
| Fixed | Original | 0.20849953 | 0.02120638 |
| Fixed | Relearned | 0.20731249 | 0.02129993 |
| Refit, grouped CV | Original | 0.20052300 | 0.02063135 |
| Refit, grouped CV | Relearned | 0.19665559 | 0.02067566 |

The archived row-CV last-row SE was 0.042989 on the same row map. That inflation
is absent in this corrected draw; this does not establish that grouping alone
explains all historical undercoverage or that the bootstrap is valid.

Independent reconstruction errors are zero for the point estimate, direct
site-centered variance, arm contrast, Wald statistics, penalties, and rerun
optimizer weights. The variance API differs by at most 1.08e-19, fold aggregation
by 5.55e-17, and normalized signed-unconstrained KKT residuals are at most
4.24e-8 (tolerance 1e-5). Observed source weights span 0.216993--0.365913 and
target anchors 0.399622--0.423102; positivity is observed, not imposed.

Evidence is under `full_refit_grouped_cv_v19_audits/` in
`seed_000001_rho_0_draw_0001/`,
`seed_000001_rho_0_draw_0001_variance_trace_v1/`, and
`seed1_rho0_draw0_draw1_optimizer_v1/`. All payload checksums passed.

## Predeclared rho=0 expansion

After the identity and first nonidentity gates passed, remaining rho=0 draw
IDs 2--20 were submitted as array `18217670`, with `--array=2-20%3`.
Tasks 2, 3 and 4 were verified RUNNING; the remaining tasks were pending at
the array throttle. Alongside the existing rho=1 first draw `18213685`, this
respects the four-active-job ceiling. Submission is not completion or an
audit pass. The rho=1 remaining draws have not yet been submitted.

The scheduling adapter now requires both prerequisite audits to explicitly
report `origin_grouped` and a positive actual partition-record count, besides
the checksum, source/draw path and v19 package gates. A stubbed-worker positive
test and rehashed wrong-protocol/zero/missing/nonnumeric-record negative tests
passed without starting any fitting job. Adapter SHA256:
`c83b140f90dd46d03fdd39e1a508333eb69ae58a2da970ecfdffecc4ff2ce53b`.
This scheduling-only file is outside the unchanged scientific worker hash.

## First corrected nonidentity draw: seed 1, rho 1, draw 1

Job `18213685` completed 0:0 in 00:46:31. The current intermediate auditor
passes with 145 explicit origin-grouped CV records, point reconstruction error
zero, variance error 1.08e-19 and raw-moment error 2.08e-17. Final nuisance fits
have zero recorded nonconvergence, line-search failures, degeneracy, support
floor events or warnings. The independently assembled arm contrasts agree
within 1.78e-15; all audit and variance-trace payload checksums pass.

| Nuisances | Common weights | Estimate | Within-draw plug-in SE |
| --- | --- | ---: | ---: |
| Fixed | Original | 0.19778379 | 0.02512135 |
| Fixed | Relearned | 0.18540400 | 0.02614141 |
| Refit, grouped CV | Original | 0.19420232 | 0.02436824 |
| Refit, grouped CV | Relearned | 0.18347724 | 0.02498284 |

The fitted source weights span 0--0.398804; target anchors span
0.597386--0.654765. These are observed values, not newly imposed constraints.
The fitted variance is 0.000624142237585151, with target/source-1/source-2
contributions approximately 0.0003650784, 0.0000008036 and 0.0002582603.
These are within-draw plug-in variances, not across-draw SDs or evidence of
sampling coverage. Evidence directories under `full_refit_grouped_cv_v19_audits/`:
`seed_000001_rho_1_draw_0001/` and
`seed_000001_rho_1_draw_0001_variance_trace_v1/`.

The independent origin audit confirms all 145 records and exact row-map and
multiplicity pairing with rho=0 draw 1. The optimizer auditor initially failed
because it hardcoded the rho=0 source reference, not because of an estimator
failure. Its reference selection now uses the actual rho and requires matching
simulation/rho/source and identity/positive-draw metadata, including CSV/RDS
agreement. Eighteen boundary/selection assertions pass. Fresh rho=0 and rho=1
regressions both pass under auditor SHA256
`f00838289f1faebe6c650f333e299359de4f08b6c0cf860e8845f408a37a0487`.
Use `seed1_rho0_draw0_draw1_optimizer_v2/` and
`seed1_rho1_draw0_draw1_optimizer_v2/`; older v1 outputs remain unchanged.
The rho=1 maximum normalized signed KKT residual is 3.68e-8, and exact zero
source weights satisfy the L1 subgradient conditions. No audit threshold was
weakened and the scientific worker/package fingerprints are unchanged.

## Predeclared rho=1 expansion and subsequent rho=0 draws

After all prerequisite reviews passed, rho=1 remaining IDs 2--20 were submitted
as array `18218871` with `--array=2-20%1`. Task 2 was verified RUNNING with
`ArrayTaskThrottle=1`; alongside rho=0 array `18217670` at throttle 3, the
two arrays have at most four active full-refit jobs. The complete predeclared
20 nonidentity draws per rho remain required; these are not new independent
simulation datasets or a passed inference gate.

Rho=0 draws 2, 3 and 4 subsequently completed 0:0 in 00:50:04, 00:44:45 and
00:35:39 respectively. Draw 4's general audit passes with 145 explicit valid
partitions, zero point/variance reconstruction error, zero final-fit failures
or warnings, estimate 0.211228716836282 and variance 0.000473557991647253.
Draws 2 and 3 subsequently passed the same auditor: each has 145 valid explicit
CV records, zero point error, variance error at most 1.08e-19, and no final
nuisance failures or warnings. Their fitted estimates are 0.208382476828303 and
0.173266606284666, with variances 0.000425942546159597 and 0.000487185786263739.
No statistical conclusion is drawn from this early subset, and pending IDs
are not replaced or omitted.

Independent paired variance traces for draws 2--4 confirm identical grouped
versus archived row-CV multiplicities and target row-map/fold/A/Y fields.
Fixed-nuisance A/B results also match exactly. The refitted plug-in SEs are:

| Draw | Grouped CV, original weights | Grouped CV, relearned weights | Row CV, original weights | Row CV, relearned weights |
| ---: | ---: | ---: | ---: | ---: |
| 2 | 0.02063388 | 0.02063838 | 0.06891024 | 0.06853161 |
| 3 | 0.02250935 | 0.02207229 | 0.03738743 | 0.03849082 |
| 4 | 0.02191117 | 0.02176139 | 0.03814136 | 0.03851743 |

The refitted target treated-arm maximum absolute AIPW values are approximately
3.137, 7.040 and 5.734 under grouped CV, versus 76.964, 16.169 and 18.240 under
row CV. All three trace reconstruction gates and payload checksums pass. These
same-resample contrasts localize the instability in the archived protocol;
they do not establish bootstrap or independent-sampling coverage.

The explicit early subset summary `rho0_summary_draws_0_4_v1/` passes all
checksums and excludes identity draw 0. Across its four positive draws, A/B/C/D
SDs are 0.01325732, 0.01433246, 0.00955706 and 0.01727061. The paired-change
variance decomposition error is zero and the full-D decomposition error is
5.42e-20, including all covariance terms. Its `diagnostic_draw_set_complete`
flag refers only to the requested subset 0:4, **not** completion of the
predeclared 20-draw pilot. `inference_validated=FALSE` and `small_B_warning=TRUE`.

## Complete rho=0 20-draw computational checkpoint

All predeclared positive rho=0 draw IDs 1--20 completed successfully and passed
the general intermediate auditor. The 2,900 actual nuisance-CV records are
explicit and valid. Maximum errors: point zero, variance 1.08e-19, raw moments
3.47e-17. Recorded final nonconvergence, line-search failure, outcome degeneracy
and support-floor counts are all zero. One warning is retained for review;
draw 13's larger plug-in variance (0.0009442926 versus roughly 0.00042--0.00049
for the other draws) is receiving a targeted intermediate-value trace.
These flags are not silently discarded because reconstruction passes.

The complete summary `rho0_summary_draws_0_20_v1/` passes its payload hashes,
excludes identity draw zero, and records 20/20 successful attempts. Across-draw
SDs for A/B/C/D are 0.01776742, 0.01859010, 0.01450047 and 0.01816599. Both
the full-D and D-minus-A covariance decomposition errors are zero. Negative
covariance terms materially offset positive component variances, so summing
only weight/nuisance variance terms is not a justified correction.

This completes one conditional resampling pilot at rho=0, not independent
Monte Carlo calibration, final variance selection, or the full experiment.
The summary retains `inference_validated=FALSE` and `small_B_warning=TRUE`.

## Remaining rho=1 batch and targeted tail/warning review

Rho=1 draws 1--15 have completed and passed their general intermediate audits.
The explicit 0:14 subset summary `rho1_summary_draws_0_14_v1/` passes hashes,
excludes identity zero and reports 14 successful positive draws. A/B/C/D SDs
are 0.02041922, 0.02568782, 0.02036135 and 0.02917242. This remains an early
subset, not the final 20-draw checkpoint or coverage validation.

With rho=0's array fully completed, the existing rho=1 array throttle was
raised from 1 to 4 without adding or replacing any draw. `scontrol update`
returned a nonzero message for already-finished task 15, but read-only checks
confirmed that both live task 16 and pending task 20 have
`ArrayTaskThrottle=4`, with tasks 16--19 RUNNING and 20 waiting at the cap.
Thus the change took effect despite that partial-command status; it was not
retried and no task was resubmitted. The total four-active-job ceiling holds.

The retained glmnet warning is in draw 12 for both rhos: error -80, convergence
not reached for the 80th lambda after 1,000,000 iterations, with larger-lambda
solutions returned. It is distinct from draw 13's high variance. The corrected
target-CV replay below locates it precisely rather than classifying it from
the final solver counters alone.

The draw-13 variance trace passes reconstruction and hashes, but confirms a
real tail: source-1 original row 776 appears twice in the resample; its
centered, site-scaled influence grows from about 5.9 under fixed nuisances to
about 42.4 after refitting. Source 1 contributes about 57% of the final variance,
and its top five positions account for about 75% of its site contribution.
The target treated-arm maximum absolute AIPW is about 15.7, with no probability
floor hit; target variance is only about 25% of the total. Actual grouped CV
is valid at all 145 recorded calls and draw 13 has no warnings. The source
influence factor decomposition subsequently reproduces the saved values and
identifies the density ratio increase (18.32 to 131.27) as the principal
amplifier, with the residual slightly smaller and weight only slightly larger.
The density logit -4.877 remains inside M_tau=5. See the explicit paired-factor
table in `STEPWISE_FUNCTION_REVIEW.md`; intact groups and a satisfied fitting
bound do not guarantee light tails or valid unconditional inference.

## Exact draw-12 target-CV warning replay

The versioned output `draw12_target_glmnet_warning_replay_v2/` passes both
payload checksums. It replays only target tuning, not source fits, aggregation
or a full refit. All 75 target calls (25 PS and 25 OR per arm) reproduce saved
prediction lengths and values exactly. Their caller, complete origin-group
vectors and actual CV fold IDs each uniquely match a saved partition record;
maximum prediction and fold-ID errors are zero.

The exploratory replay initially hardcoded Gaussian OR fits, whereas the
saved reference family is Binomial. That diagnostic-script defect was corrected
to use `reference$family` and its resolved GLM family. The old CSV in the
diagnostic source directory is retained but is **not** formal replay evidence.
The corrected script fingerprint is
`e8ccf02ea660f6477ec5ab0aba73e9c07f319dd6c0e06b8b50202ca22252ca82`.

Exactly one of the 75 corrected calls reproduces the warning: target PS,
outer fold k1=5, k2=NULL, 801 training rows and 199 held-out prediction rows,
five origin-grouped CV folds. The returned path length is 75;
lambda.min=0.0495189105881761 and lambda.1se=0.0718434605610932. The specified
minimum-loss rule selects lambda.min and its 199 propensity predictions match
the saved fit exactly. The other 74 calls have no warning. Both rhos share this
target resample and seed, so no duplicate full refit was needed.

This establishes the warning's source and reproducible finite predictions;
it does not show that the unavailable lower-penalty path could never change
the tuning choice. The saved replay does not include the full lambda vector
or internal `jerr`, so neither a selected path index nor an exact internal
truncation relationship is asserted. The warning remains part of the record.

## Final paired conditional-pilot comparison

The complete summaries `rho0_summary_draws_0_20_v1/` and
`rho1_summary_draws_0_20_v1/` both pass their payload hashes, include exactly
20 positive draws, and exclude identity zero. A/B/C/D denote fixed nuisances
with original/relearned weights and refitted nuisances with original/relearned
weights, respectively.

| Rho | SD A | SD B | SD C | SD D | Original analytic SE |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.01776742 | 0.01859010 | 0.01450047 | 0.01816599 | 0.02169048 |
| 1 | 0.02110382 | 0.02778184 | 0.02021331 | 0.02823550 | 0.02539934 |

For rho=1, Var(D)=0.000797243234690575; both the D-minus-A and complete-D
signed covariance decompositions reconstruct with zero error. In particular,
twice the covariance between baseline A and its total change is approximately
+0.0000878504 at rho=1, versus -0.000135078 at rho=0. Its direction is not fixed,
so separately adding positive variance components is not a valid correction.

D-minus-original-reference mean shifts are -0.007345 at rho=0 and +0.008159 at
rho=1, with conditional Monte Carlo SEs about 0.004062 and 0.006314. Twenty draws
give noisy SD ratios and do not establish bias, coverage, or method superiority.
The full-refit pipeline is now computationally auditable; the primary analytic
CI remains unchanged. Further repetition of this one dataset would improve
conditional Monte Carlo precision but still would not answer independent-
sampling coverage. The next stage therefore uses a predeclared fresh-seed
inference pilot, not automatic promotion of either conditional SD to paper SEs.
