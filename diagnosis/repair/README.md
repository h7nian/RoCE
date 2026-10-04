# Repair implementation in progress

The user authorized implementation of the reviewed route, including fold-summed
source outer/inner and target outer/inner calibration; pooling all inner-fold
moments before learning each outer fold's weights; one-/two-round message
contracts; and selectable common-TATE, separate-arm and joint-TATE aggregation.
The new nested validation has three fold roles. Communication rounds remain a
separate property to verify against the complete inference message contract.

Following the user's 2026-09-20 clarification after the interim review,
performance assessment prioritizes `joint_tate`, compares it with
`separate_arms`, and retains `common_tate` as a constrained comparison.
High-dimensional expansion does not require every aggregation option to
perform well. Common-weight coverage or RMSE alone does not block expansion;
shared implementation defects and inference validity still require review.
All variants, unavailable repeats and Monte Carlo uncertainty remain reported.

Results and build products belong under
`/scratch.global/zhan9381/FACE-HD/implementation/`. Code stays in this repository.
Do not install over any frozen v4 library.

## Current R0 changes

- Initial nuisance fits use deterministic training-subset seeds and exact-input
  caching. No cross-subset lambda or initial warm-start sharing.
- Final nuisance tuning uses its own deterministic training-subset stream.
- Target prediction/full-fit caches check actual data before reuse.
- Rho reuse records and validates cache/training policy and solver metadata.
- Final C++ outcome fitting uses the common Newton/CD dispatcher.
- `nuisance_solver` arguments select scoped C++ solver overrides; both solvers
  remain available. glmnet fits are labeled separately.
- Unsupported cross-fitted `use_rcal=TRUE`, nonfinite supplied penalties and
  dimension-error fallbacks no longer silently change behavior.

The shared loss-assembly engine now serves source and target calibration.
`target_nuisance_method="hou_calibrated"` adds calibrated outer and inner
anchors. `source_validation_method="calibrated"` adds complete inner source
training with three excluded fold roles. Defaults preserve the legacy options
for explicit component comparisons. Both calibrated inner options require at
least four folds; the study design remains ten folds and 1000 observations/site.

`aggregation_mode` selects `common_tate`, `separate_arms`, or `joint_tate` in
the TATE API and reaggregation API. Joint optimization uses full within-site
two-arm covariance, the actual number of observations, and arm-specific
discrepancies. Its empirical-mass derivative includes both coordinates at each
source. This is a conditional calculation, not a proof of full nuisance or
selection inference. The older multiplier-bootstrap routine explicitly rejects
two-vector results until its reconstruction has been extended.

The simulation APIs propagate solver, nuisance method and matrix layout.
`additional_aggregation_modes` adds paired alternative rows from the same
nuisance fits. The original ordinary target-only TATE remains the fixed
benchmark when a calibrated RoCE anchor is selected.

## Validation state

- `implementation/r0/before_source/` preserves the pre-edit package source.
- `implementation/r0/r_only_checks/` passed three groups / 49 expectations for
  R caching, RNG, explicit errors and initial-fit isolation (both arms and both
  protocols, a 1000/site small-basis fixture). This used the frozen native
  library and **does not validate the new C++ dispatch or solver override**.
- Slurm job **1428325** was canceled while pending: its source stage was
  superseded. It produced no validation result.
- `implementation/r0/checks_v2/` built the native code on the current node.
  Selected solver, loss-assembly, caching, CV and model tests passed except
  one old parallel-equivalence fixture with only 90 observations per site:
  its independently tuned initial density-ratio CV had no converged lambda.
  The failure log is retained; grid feasibility is not thereby resolved.
- `implementation/r0/rho_equivalence_1000/` reran that equivalence group with
  1000 observations per site and passed all 67 expectations, including actual
  PSOCK execution. The v2 native tests also verified explicit Newton/CD
  selection, independent KKT checks, and compact/block calibration equivalence
  including selected penalties. These are component checks, not MC coverage
  validation or a complete package check.
- Small checks may run on the current node, as explicitly authorized. Always
  compare the stage manifest to the source actually compiled.
- `implementation/r1/checks_v1/` passed native compilation and all selected
  calibration, nested-isolation, solver, cache and rho-reuse groups. The new
  1000/site adapter tests independently reconstruct both target loss gradients
  and perturb both excluded folds under both protocols and arms.
- `implementation/r1/checks_v2/` built the first three-mode API. Integration
  and derivative tests ran, but the suite reported six test assertions: one
  name-attribute mismatch, four uses of an unavailable named source index,
  and one overly strict relative comparison of a small finite difference.
  Independent printing found derivative absolute differences below 6e-10;
  corrected tests now require absolute error below 1e-8. The failed stage is
  retained, and subsequent changes require a fresh stage.
- `implementation/r1/checks_v3/` passed the corrected joint derivative,
  three-mode reaggregation, complete inner-calibration, and simulation comparison
  tests, including preservation of the fixed target-only TATE reference.
- `implementation/r1/full_suite_v3/` ran 348 test cases: one obsolete
  source-comment assertion and four erroring tests prevented a full pass.
  One error was an aggregation KKT mismatch; three arose in 120/site fixtures
  with no admissible initial-CV penalty. The source-comment test was replaced
  by the existing independent objective tests. Weight-derivative fixtures now
  use the requested 1000/site scale; the original numerical failure is retained
  separately and has its own unbounded-objective regression test.
- `implementation/r1/initial_cv_failure_7102/` reproduces the original failure.
  In CV fold 2, the original largest penalty 0.48946 is insufficient: an
  independently solved recession direction has nonnegative source projections
  and slope -0.04192, proving no finite minimizer for that penalty or any
  smaller one. The old CD reported convergence near coefficient magnitude 100
  with KKT 0.0551, while Newton rejected the fit. A larger null-model penalty
  yields a stationary intercept-only fit in both solvers. No production lambda
  grid has been silently enlarged. Whether p100 study cells need an explicit
  grid-boundary repair remains to be checked.
- The R1 repair certifies native nuisance convergence by KKT, applies
  the same coefficient-bound rejection to both solvers, and checks the common
  aggregation optimizer's KKT before accepting convergence. Nuisance failures
  abort the calibrated training task with site/fold context. Training metadata
  became `training_subset_v2_kkt`; the following stages tested this change.
- `implementation/r1/checks_v4/` built and passed the new nuisance KKT tests,
  including rejection of a coefficient-bound pseudo-solution, but retained one
  aggregation error. That residual came from separately flooring variance
  components in the optimizer while differentiating their unfloored values.
  The optimizer now preserves empirical variance components and floors only
  the Wald discrepancy variance, matching the derivative's objective.
- `implementation/r1/checks_v5/` passed the selected regression groups after
  this correction, including the formerly failing untruncated-inference case.
- `implementation/r1/acceleration_v4/` compared CD/no cache/block/serial with
  Newton, caching, compact matrices and parallel execution at 1000/site,
  p4, four folds, K2, binary outcomes. Times were 23.9/8.62/2.63/2.64/1.78s;
  the largest TATE difference was 4.87e-8 and largest SE difference 4.63e-9.
  The strict overall gate failed because upstream solver differences moved
  some data-dependent lambda values beyond the prespecified 1e-10 tolerance.
  This is not labeled a complete equivalence pass.
- `nuisance_tol` now explicitly controls final-fit precision. The first
  `acceleration_v5` rerun at 1e-10 correctly stopped on a Newton KKT residual
  around 1e-9. Near stationarity, subtracting two complete L1 norms erased a
  small predicted decrease. The working tree now computes that decrease by
  combining coordinatewise smooth and L1 slopes before multiplying by the
  step; a focused regression test and new native validation are pending.
- `implementation/r1/checks_v6/` passed those focused checks.
  `implementation/r1/full_suite_v6/` then passed 351 test cases / 2770
  expectations, with no errors or failed expectations (120 recorded warnings,
  four skips). Package help was regenerated in a separate scratch tree.
- `implementation/r1/acceleration_v6/` passed the original strict comparison
  gates at explicit `nuisance_tol=1e-10`: all compared penalties agreed within
  1e-10, maximum TATE difference 4.90e-11, maximum SE difference 5.67e-12,
  maximum weight difference 7.57e-10. Times were 24.56/8.70/2.64/2.64/1.81s.
  This remains a p4/four-fold numerical fixture, not a coverage study.
- `implementation/r1/large_basis_v6/` passed a deepest-layer initial weight
  fit at raw p100 / 200 working features, 1000/site and ten original folds.
  The selected source training set had 697 rows / 275 arm observations.
  CD and Newton selected the same lambda 0.1774232701; maximum prediction
  difference was 4.90e-10. Times were 144.5s and 4.37s. Both rejected 87 of
  100 CV candidates. This validates one high-dimensional nuisance task only.

An independent statistical check is now in
`diagnosis/repair_audit/clipped_score_review.md` on scratch, with eight
population rows under `clipped_score_checks_v1/`. The current clipped plug-in
losses can leave nonzero nuisance sensitivities even with one model correct.
Matching the OR derivative and tilt derivative to the actual evaluation score
eliminates these sensitivities in the independent examples. The corresponding
recipe is now available through the experimental `calibration_control` argument:
`recipe="score_derivative"`, `target_propensity_initialization="calibrated"`,
and explicit `target_radius=log(9)` for the current clipped target propensity
DGP. These change the statistical procedure; they are not acceleration tricks.
Defaults retain the legacy recipe. Higher-order/rate, selection, baseline and
communication arguments remain open. A numerical test pass is
not a declaration that the full statistical repair is complete.

The v7 focused checks passed exact final-calibration fit reuse and simulation
checkpoint/configuration and CPU-budget guards. Its full package check had one
NOTE for an unqualified `plogis`; the working tree now uses `stats::plogis`.
The p100 v7 component check completed both full inner source/target fits;
only its postprocessing needed a namespace correction, recorded in the output.

`implementation/r2/checks_v2/` built the new native recipe and passed 359 test
cases / 2892 expectations, with 120 warnings and four skips. This includes
independent score-gradient checks in four model-correctness branches and a
small-basis three-level pipeline. A subsequent correction classifies the two
new aggregation rows as TATE in downstream summaries; `r2/checks_v3/` passed
the full 359 cases / 2894 expectations with the same warnings/skips.

`r2/calibration_cost_v2/` passed the score-derivative recipe's complete deepest
inner task at p100, 1000/site, ten folds, and 100 lambdas: source 50.53s, target
48.16s, plus target messages 2.06s. This checks one treated-arm task only.
`r2/acceleration_v2/` passed all original strict gates for the new recipe's
small-basis three-level pipeline: maximum TATE difference 8.91e-12, SE
difference 1.12e-12 and weight difference 3.63e-10. Selected penalties passed
the 1e-10 gate. CD/no-cache/serial took 37.54s and the accelerated path 1.73s.

The full `r2/package_check_v3/` completed its tests with one documentation
WARNING (two new arguments lacked prose). Their roxygen/Rd entries are now
fixed; `package_check_v4/` returned `Status: OK` with `--no-tests`. The full
tests had already passed for the same implementation in v3.

`stage_source.py` creates a new source snapshot. `run_checks.sh` compiles that
snapshot and runs its tests using an isolated installed library; an optional
second argument is the testthat file filter. Neither script writes build
artifacts to the working tree. The generated Rcpp registration files must be
regenerated in each stage because the scoped solver setter adds a binding.

## Shared caches and the grid-cost review (R4–R5)

`nuisance_cache_dir` optionally enables a private per-run cache shared by source
workers across outer folds. Keys include the exact fitting inputs, solver and
training policy; atomic entries verify their model digest. The owner process
cleans up its directory. `NULL` retains the memory-only cache. Supply an existing
writable directory on scratch and retain `use_lambda_cache=TRUE`.

R4 `checks_v2` passed 367 cases / 3083 expectations, with zero failures/errors,
133 recorded warnings and four skips. `package_check_v2` returned `Status: OK`
with `--no-tests`; the full regression suite was separately run on the same
implementation. Strict small-basis acceleration gates passed: maximum TATE
change 8.91e-12, SE change 1.12e-12 and weight change 3.63e-10. Shared disk caching
was slower on that small fixture (8.45s versus 1.63s with memory caching).

For one complete deepest source calibration task at p100 / 200 features,
1000/site, ten folds and 100 lambdas, separate worker waves took 51.90/52.15s
with memory caches and 53.03/0.93s with a shared cache. The final coefficients
and selected penalties were bit-identical. This component reuse result does
not establish the full pipeline speedup. Full C3/K2/rho0/seed1 job 1447194 uses
R4, shared caching and the score-derivative recipe. Its reference is the same
recipe in R2 job 1437467; their paired p100 check remains pending.

The pilot driver now divides its total allocated CPU budget between concurrent
arms. Every repeat remains a separate job, with a single multi-partition
request. The corrected driver completed both recipes on the existing small-basis
smoke fixture and preserved the ordinary target-only reference exactly.

R5 changes diagnostic scripts only. Its p100 component study finds that native
CV dominates final-refit time, but requested grid size is not a universal cost
predictor: initial source-weight CV took 3.90/3.53/6.24s for 100/50/25 lambdas.
Of 500 nominal fold/penalty combinations at 100 points, 416 were already
skipped; 15 failed attempts accounted for about 96% of measured fit time.
Calibrated OR does benefit from fewer grid points, but selected penalties and
predictions change. Reducing both final source calibration grids to 50 changed
one held-out mu1 by -0.01254 with fixed initial models. This is not TATE or
coverage evidence, and the running/default 100-point configuration is retained.
See [grid review](/scratch.global/zhan9381/FACE-HD/implementation/r5/lambda_grid_review.md).

`compare_repeat_pilots.R` checks frozen pilot pairs, including input identity,
training subsets, all saved nuisance penalties/coefficients, TATE rows, and
weights from all three aggregation modes. Missing/failed recipes yield an
explicit incomplete result. The completed R2/R4 small-basis pair passed before
scheduling the p100 review; altered-model/subset guard checks are separate.

The subsequent live p100 cache inspection found slow **control-arm final OR
CV** tasks: source s2 outer fold1 took 785.83s and target outer fold1 took
289.82s, both with successful selected fits. Their exact initial/final fitting
inputs were reconstructed and verified against every cache hash without
refitting. `capture_cached_calibration.R` and `check_outcome_grid.R` preserve
and profile these cases. The profiler's eight-penalty prefix matches native
CV scores within 1e-12; full-path diagnostic job 1449410 is queued. See
[slow OR evidence](/scratch.global/zhan9381/FACE-HD/implementation/r5/slow_or_review.md).
The initial-weight timing fraction must not be extrapolated to the full method.

Dependent jobs 1448996 and 1448997 review the accelerated task and full p100
pair respectively. They run after termination, including failure/timeout;
missing recipes are explicitly incomplete. The
[validation handoff](/scratch.global/zhan9381/FACE-HD/implementation/r5/validation_status.md)
records jobs, output paths, numerical gates and the next expansion criteria.

## Remaining implementation priorities

1. Finish the full p100 shared-cache/reference pairing, retaining numerical
   gates and failed fits; then expand difficult cells with fixed seeds 1–5–10.
2. Add more granular task timing and review failed CV attempts and grid bounds.
   Keep grid experiments distinct from result-preserving acceleration.
3. Complete clipping/nuisance/selection inference, message-only reconstruction
   and unequal-fold scaling. Baseline conditional derivatives are checked;
   learned IVW weights and selection transitions remain open.
4. Complete the full function review and scenario comparisons before formal
   fixed-n MC500 validation.

Acceleration must preserve the statistical procedure: compare serial/no-cache/
CD with Newton, exact caching, parallel execution and compact matrix layouts on
identical data, folds and tuning grids. Check selected penalties, objective/KKT,
point estimates, SEs and weights. Statistical repairs may change results; their
effects must be reported separately from pure acceleration differences.

## Single-repeat Slurm pilot

`submit_repeat_pilot.py` prepares a fixed 1000/site, p100, ten-fold development
manifest. Each job handles exactly one scenario and seed. Within it, requested
calibration recipes share the same generated data; each recipe compares the
three aggregation modes using its own shared nuisance fits. The ordinary
target-only TATE is checked for exact preservation across recipes.

```bash
python3 diagnosis/repair/submit_repeat_pilot.py prepare \
  /scratch.global/zhan9381/FACE-HD/implementation/r2/pilot_next \
  --check-root /scratch.global/zhan9381/FACE-HD/implementation/r2/checks_v3
python3 diagnosis/repair/submit_repeat_pilot.py submit \
  /scratch.global/zhan9381/FACE-HD/implementation/r2/pilot_next \
  --max-jobs 3 --test-only
python3 diagnosis/repair/submit_repeat_pilot.py submit \
  /scratch.global/zhan9381/FACE-HD/implementation/r2/pilot_next --max-jobs 3
```

The default partition list is `saffo-2tb,msismall,agsmall,amdsmall`; Slurm selects
one partition for each job. Partition/account/resource settings are arguments.
The manifest, configuration, installed library and copied workflow are hashed.
Submission receipts and an exclusive task directory prevent automatic duplicate
submissions or concurrent execution of a repeat. Failed/ambiguous submissions
remain recorded and require inspection before recovery. All logs, models and
temporary files stay under the pilot directory on scratch. This bounded pilot
does not include unrepaired Federated-DR/Pooled-DR baselines and is not MC500.

The current pilot at `implementation/r2/pilot_v2/` contains C1/C2/C3, K2, rho0,
seed1. Jobs **1437464**, **1437466**, **1437467** were submitted on 2026-09-20
using the multi-partition list above. These three repeats are pipeline checks;
they cannot estimate scenario coverage reliably. Scheduler receipts and
task-local status files are the current execution record.

Dependent review job **1440338** runs after all three current pilot jobs end,
including failures. It writes task status, TATE estimates and paired recipe
differences to `implementation/r2/pilot_v2_review/results/`. Both successful and
failed operational fixtures were used to check this summary script. It does
not turn three development repeats into a coverage/RMSE study.

The first submission (v1 jobs 1436882-1436884) was canceled while pending after
the operational driver check found that `n_bootstrap=0` violates the shared API.
The corrected driver retains 5000 for this unused comparison-method setting.
`driver_smoke_v4/` then completed both calibration recipes and all three
aggregation modes at 1000/site with the smaller operational-test basis/folds,
including exact target-reference equality and artifact saving. The failed
`driver_smoke_v3/` record is retained; neither operational check is a p100 study.

## R3 baseline inference repair

The new `R/baseline_influence.R` differentiates the actual fitted comparison
models conditional on their selected penalties and active sets. MLE nuisance
contributions use the correct positive sign; lasso contributions include the
empirical feature-scale derivative. Raw predictions enter likelihood scores,
while clipping indicators enter the evaluation sensitivities. The constant OR
fallback is differentiated as its empirical mean. Retained metadata does not
change the ordinary target-only predictions or RNG stream.

Density corrections use the raw fitted tilt, account for empirical normalization
and clipping, and retain active-coordinate metadata. Pooled-DR computes its own
sensitivities at the pooled prediction/estimate and uses the correct source and
target sample-size factors. Its unused source-local AIPW fits and
`dr_site_components` argument were removed. Pooled nuisance fitting now has its
own arm-specific seed, recorded in the result; historical pooled penalties and
point estimates can change. The change is versioned as
`baseline_inference_policy="conditional_active_set_v1"` in simulations and
checkpoints. This is a baseline repair, not an acceleration-equivalence claim.

`r3/baseline_v1` passed 45 targeted cases / 375 expectations; `baseline_v2`
passed 47 / 422, including complete case-weight refits for Pooled-DR. The full
`checks_v3` suite's first run passed 364 cases / 3060 expectations, with 133
warnings and four skips, but its editable workspace shell wrapper was changed
while active and exited after tests. Those results are archived separately;
the immutable staged driver rerun completed normally with the same 364 cases /
3060 expectations passing. Always invoke the frozen script
at `CHECK_ROOT/source/diagnosis/repair/run_checks.sh`, not a workspace script
that may be edited during execution. New drivers also verify the stage manifest
before compilation.

`r3/package_check_v3` returned `Status: OK` with tests separately checked.
`r3/p100_baselines_v3` completed a C3/K2, 1000/site, p100 fixture in 103.5s and
verified paired-arm variance reconstruction. This does not establish coverage.
See the [baseline derivation](/scratch.global/zhan9381/FACE-HD/implementation/r3/baseline_inference.md)
and [validation record](/scratch.global/zhan9381/FACE-HD/implementation/r3/validation_status.md).
The fitted IVW-weight derivative, selection transitions, transport assumptions,
and high-dimensional remainders still gate formal baseline inference claims.
