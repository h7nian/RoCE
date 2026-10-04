# Target calibration and production audit (2026-09-19)

The evaluation design is fixed at **1000 observations at every site**, raw
dimension 100, the existing 200-column `[X - kappa, X^2]` working basis, and ten
outer folds. Larger-sample experiments are diagnostic evidence, not a remedy.
The primary endpoint is the direct, common-weight TATE; armwise TATE is reported
separately. Nominal coverage is 95%. All Monte Carlo comparisons retain counts,
uncertainty, failures, and paired comparisons against the target-only anchor.

Executable diagnostic code lives here. All results, installed test libraries,
package-check stages, and saved fits live at:

`/scratch.global/zhan9381/FACE-HD/diagnosis/fixed_n_audit_20260919/`

The production library is read-only throughout this diagnostic:
`results/direct_tate_mc500_b5000/Rlib_production_20260917_v4`.
Nothing here changes `R/`, `src/`, the production manifests, or production raw
results. Adding a diagnostic alternative is not adopting a new production method.

## Experiments and their scope

`run_diagnostic.R` has three explicit modes:

- `smoke`: first evaluation fold of each arm; checks implementation and runtime.
  Its numerical results are **not** full-target estimates or performance evidence.
- `target`: all outer target folds, using the real production complement-fold
  estimator, the actual FACE DGP, and identical data/folds for each alternative.
- `full`: first fits the complete one-round RoCE estimator. Replaces the target
  anchor in **both outer and inner folds**, then reruns the existing arm and TATE
  aggregation. Source nuisance fits and one-round initial outcome summaries stay
  fixed. A baseline reconstruction must reproduce estimates, SEs and weights.
  This tests changing the target anchor, not changing every target-side fit.

Every mode compares six variants: production lasso; oracle outcome; oracle
propensity; double oracle; calibrated propensity with the original outcome fit;
and calibrated propensity with its matched weighted outcome fit. For each held-out
target observation it also records the conditional product remainder
`(1 - p_true / p_fitted) * (m_fitted - m_true)`. This identifies the observed
target remainder without treating a big-O bound as a measured error.
In particular, `oracle_both` means **both target-anchor nuisances** are replaced;
it never means the entire RoCE estimator has oracle source nuisances. New driver
outputs also label this explicitly with `intervention_scope=target_anchor_only`.

Two RCAL backends are available:

- `reference` calls installed RCAL 2.0 directly. Features are standardized on
  complement-training observations; the outcome fit uses arm-specific odds
  weights `(1 - p_a) / p_a`. Both arms get separate calibrated propensity fits.
- `native` reuses RoCE's existing compiled exponential-tilting solver for the
  identical calibrated propensity objective, and glmnet for the matched weighted
  logistic objective. No new optimizer is implemented. Fixed-penalty solutions
  are checked against RCAL; selected fits must pass a KKT check. The two backends
  need not pick identical penalties because their CV partition/path implementations
  differ. They are not described as bitwise-equivalent complete estimators.

The declared propensity grid has 100 candidates over a relative range of 1e-4.
CV fitting allows 100 iterations, final fitting 1000, with convergence/KKT checks.
Invalid candidates are excluded, with path-tail stopping at the first failed fit;
counts are saved. The reference package uses a similar tail exclusion policy.
Native weighted glmnet uses up to 100 candidates and its normal path stopping.
Its per-candidate nonconvergence vector is unavailable and is reported as NA;
warnings are recorded separately. Trial 1 originally wrote a placeholder zero
in that diagnostic field; it must not be interpreted as verified convergence
of every outcome-CV candidate. This metadata-only correction does not change
fits, estimates or SEs. The exact executed scripts are retained in scratch at
`code_snapshot_trial1/` with SHA-256 hashes.
Propensity clipping matches the existing target estimator. Training scaling,
fitted probabilities, convergence, warnings and code hashes are recorded.
Outer folds are shared across candidates. Nuisance-CV partitions follow each
backend's documented seed/partition recipe; they are not asserted identical
across glmnet, the RCAL reference package and the native calibration backend.

These are **candidate methods**. Neither the package name RCAL nor an improvement
in a small pilot establishes nominal coverage in C3, where the propensity working
model is misspecified and the outcome is logistic rather than linear.

## Reproduction

Load the cluster R module and set `R_LIBS_USER=/users/0/zhan9381/Rlibs`,
`ROCE_PROJECT_LIB` to the frozen library, and `ROCE_DIAGNOSTIC_ROOT` to the
scratch directory. `run_diagnostic.sh` requires the project root, mode and
configuration explicitly; the Slurm array index is the seed. Existing seed
directories are never overwritten. Only directories containing `COMPLETED`
are successful replicates; absent completions must not be silently dropped.

Validation scripts:

- `audit_functions.R`: parse all package, launcher, workflow and test R files;
  inventory top-level functions, parameters, package calls and possible unused
  parameters. Text mentions in tests do **not** mean behavioral test coverage.
- `test_target_rcal.R`: AIPW identities, treatment-label symmetry, evaluation-data
  exclusion, reproducibility, and the legacy ignored RCAL flag.
- `test_native_rcal.R`: fixed-penalty objective equivalence with RCAL.
- `test_face_protocols.R`: FACE C1/C2/C3 integration checks for both protocols
  and exact baseline reaggregation, at 1000/site with a small feature basis.
- `check_baseline_influence.R`: numerical case-weight derivatives of independently
  refitted MLE nuisance models, checking the baseline IF correction signs.
- `check_inner_fold_tuning.R`: perturb only an inner evaluation fold's outcomes
  and compare the nuisance coefficients trained on the unchanged other fold,
  with and without lambda sharing, in both communication modes.
- `check_solver_dispatch.R`: compare compiled fitting behavior in separate fresh
  R processes. The solver setting is cached on first use, so changing its
  environment variable in an already initialized process is not a valid comparison.
- `check_calibrated_kkt.R`: independently reconstruct the source calibrated-loss
  gradients from saved coefficients and data at p=100, without using the package's
  loss or gradient routines.
- `check_dgp_structure.R`: distinguish outcome/PS signal sparsity from the
  all-coordinate covariate shift and inspect fitted nuisance support.
- `analyze_full_target_tails.R`: inspect held-out inverse-probability weights in
  the saved full-pipeline diagnostic.
- `summarize_production.py`: independently recompute current v4 TATE cell metrics
  and paired MSE differences, verifying n/site, dimension, folds, finite outputs,
  coverage arithmetic, duplicate IDs and the recorded package fingerprint.
  It does not replace the production checkpoint/checksum gates. The coverage
  after Monte Carlo recentering is a diagnostic only, never a usable confidence
  interval procedure.

## Cross-fitting issue reproduced independently

The production `use_lambda_cache=TRUE` shares the first inner fold's selected
penalty with later inner folds inside the same outer fold. Although coefficient
training excludes both k1 and k2, the **selection of that reused penalty does not**:
the first tuning sample contains later inner evaluation folds.

In a FACE C1 fixture with 1000/site, p=4 and three outer folds, keep outer k1=1
and coefficient-training fold 2 unchanged; shuffle only evaluation fold 3's Y.
The recorded `per_k2_alpha$k2_3` changes as follows:

| Protocol | Cache | Lambda before | Lambda after | Maximum coefficient change |
|---|---|---:|---:|---:|
| One round, target initial model | TRUE | 0.00740768 | 0.02564255 | 0.422023 |
| One round, target initial model | FALSE | 0.00672552 | 0.00672552 | 0 |
| Two round, source initial model | TRUE | 0.01596467 | 0.02084021 | 0.130627 |
| Two round, source initial model | FALSE | 0.00971305 | 0.00971305 | 0 |

This falsifies strict **inner** training/evaluation separation for cached tuning
in both protocols. It does not show outer k1 leakage: the cache resets per k1.
Nor does this small-feature counterexample quantify its contribution to the
p=100 coverage deficit. That requires a separate paired estimator experiment.
The source initial-density-ratio loop has the same tuning-sharing pattern by
code inspection; the perturbation experiment above specifically tests initial
outcome fits. The manuscript claim and the comments that these are fully
out-of-two-fold nuisances need this qualification.

A repair should only reuse a tuning result for the **same excluded folds,
training data, arm, model and tuning settings**, with a deterministic training-set
seed; use the existing uncached path as the reference. The rho-reuse helper also
hard-codes `use_lambda_cache=TRUE`, so turning the public flag off is not enough
for a production rho-group run. Propagate and validate that setting in stored
fits and reuse guards. Do not loosen the guard or mix altered fits into v4.

## Initial confirmed findings

1. All existing regression tests pass; the isolated `R CMD check` has Status: OK.
   There are 127 test warnings, mainly small-class binomial fixtures. Passing
   these checks does not validate Monte Carlo coverage.
2. The static inventory covers 286 top-level R functions (56 exported), 111 parsed
   files and eight unused-parameter candidates. Some unused formals are explicitly
   retained for API symmetry; none is automatically declared dead or deleted.
3. `fit_site_aipw(use_rcal=TRUE, use_crossfit=TRUE)` produces the same result as
   FALSE with the same seed: that switch is ignored in the default branch.
   The old non-cross-fitted RCAL branch also omits the matched weighted outcome
   fit and uses a complemented treated propensity for the control arm.
4. The baseline IF correction sign is independently wrong for the MLE derivation:
   maximum numerical-derivative errors are 0.480 and 0.221 for the implemented
   signs, versus 2.3e-6 and 5.1e-7 for positive corrections. This does not validate
   applying either MLE formula to lasso fits; a sign-only change is not a complete
   solution to the baseline inference problem.
5. The final outcome refit still calls the coordinate-descent implementation in
   `fit_general_glm_cpp`, despite the global Newton label. This affects both
   communication modes and contradicts the old all-nuisance-Newton description.
6. Current v4 raw data include 21,654 files and 72 primary TATE scenario cells;
   counts vary from 49 to 500. Thirty-five pointwise Wilson intervals for coverage
   have upper endpoint below 0.95. Ten cells have larger RMSE than target-only;
   five have a positive paired MSE difference exceeding 1.96 MCSE. These are
   descriptive pointwise comparisons, not simultaneous tests or final MC500 gates.
7. Shared shift is a distinct aggregation problem worth isolating: C1/K4/rho=2.5
   has 500 replicates, direct TATE coverage 0.788 and RMSE 0.03484, versus armwise
   0.934 and 0.02192, and target-only 0.944 and 0.03180. A target-only intervention
   cannot be assumed to resolve this difference.
8. In that shared-shift cell, source s1 has mean screening discrepancy -0.14574
   and evaluation discrepancy -0.14268. Its mean Wald statistic is 3.195;
   penalties activate in 98.94% of folds, yet mean direct-TATE weight is 0.11551.
   Armwise treated/control weights are 0.000136/0. The observed issue is residual
   borrowing after detection, consistent with the documented soft penalty, not
   an established bug in the screening statistic. The corresponding treated-only
   shift has discrepancy +0.29125 and mean direct weight 0.00269.
9. After subtracting the Monte Carlo mean bias solely as a diagnostic,
   C3/rho=0 direct-TATE coverage becomes 0.962/0.970/0.963 for K=2/4/8;
   shared-shift C1/K4/rho=2.5 becomes 0.952. This supports centering error as a
   major contributor, without proving all variance calculations correct.
10. The armwise diagnostic currently uses joint fixed-weight pseudo-value
    variance, whereas direct TATE includes its weight-learning contribution.
    Their coverage values therefore compare complete reported procedures, not
    a pure weight-rule intervention with identical variance definitions. Their
    RMSE comparison is unaffected by this distinction. Armwise is a useful
    diagnostic, not an automatically validated replacement primary method.

## Completed RCAL target-anchor pilot

All three configurations have all 40 prespecified seeds, including the C1 seed-23
retry with unchanged estimator settings. The initial 30-minute timeout is retained
under `failed_attempts/`; a smaller-resource retry completed in about ten minutes.
These are target-anchor TATE metrics, **not full-method MC40 coverage results**.

| Config | Original bias | RCAL + weighted OR bias | Original RMSE | RCAL + weighted OR RMSE |
|---|---:|---:|---:|---:|
| C1 | -0.002302 | +0.000134 | 0.033326 | 0.034473 |
| C2 | -0.005572 | -0.002992 | 0.035688 | 0.038545 |
| C3 | -0.007560 | -0.007585 | 0.028765 | 0.032150 |

The C2 paired MSE increase is 0.000212 (MCSE 0.0000852). The C3 paired estimate
shift is approximately -0.00002 (MCSE 0.00177); its measured held-out product remainder is
-0.01181 versus -0.01103 for the original fit. This configuration provides no
evidence of a general bias/RMSE repair. MC40 is too small to certify nominal
coverage. The comparison concerns this candidate's explicitly specified tuning
recipe, not every possible RCAL implementation or tuning choice.

The full C3/K2/seed-1 integration run passed exact baseline reconstruction for
point estimate, SE and weights. Its outer target results agree exactly with the
independent target-only diagnostic. It rebuilds both inner and outer target inputs
while retaining source nuisance fits; one seed is an integration check, not a
performance comparison. In that seed, a treated observation has true propensity
0.1 but RCAL predicts 0.01031, producing inverse weight 96.99 (original fit: 17.31).
That is an observed route to increased variability.

Independent KKT reconstruction on this original full fit covers 80 final source
nuisance fits (two arms, ten folds, two sources, two models). The largest absolute
residual is 5.33e-7; the largest residual relative to its penalty is 1.80e-5.
Thus these selected fits solve their stated calibrated objectives accurately.
This does not validate every CV candidate or every production replicate.

Fresh-process solver probes also confirm the dispatch inconsistency: final
outcome fitting takes seven iterations with identical coefficients under both
requested solvers; the initial tilting positive control changes from nine to six
iterations and agrees in coefficients within 2.82e-9. Source inspection identifies
the final outcome's hard-coded coordinate-descent loop.

## Why the simple sparsity calculation is inadequate

`FACE_SIGNAL_COORDINATES=4` describes outcome and propensity signals. The source
covariate generator shifts **all 100 coordinates**. For a source skewness nu,
the exact log density ratio is a sum of `log(2*Phi(nu*(X_l-kappa)))` over all p
coordinates; the ideal transported weight therefore depends on the 96 coordinates
absent from the outcome/PS signals as well. Outcome sparsity cannot be reused as
the sparsity of the merged site/treatment nuisance.

In the saved C3/K2 fit, final tilting models average 27.7 active coefficients for
the treated arm and 67.1 for control, predominantly in outcome-irrelevant
coordinates. The log-density-ratio SD from those 96 coordinates is 0.803 and 1.621
for sources with nu=0.1 and 0.2, respectively. These observed supports do not
identify the exact population projection's sparsity, but reinforce the need to
measure its complexity rather than set both nuisance sparsities to four or eight.

There is also a literal theory/design qualification: the manuscript assumes
sub-Gaussian working features, while the working basis contains squares of
unbounded normal/skew-normal variables. Those squared coordinates are not
sub-Gaussian. None of these points proves that n=1000 cannot work; they prevent
interpreting a plug-in asymptotic rate calculation as a validated finite-sample
diagnosis.

The production holes are not all reuse-guard cases. C2/K2, K4 and K8 are each
missing seed 26 within seeds 1-50. Shared-shift C1/K2 seed 108 timed out after
12 hours (job 1241257_108). Negative-transfer C1/K8 has ten failed tasks in its
n=300 rung (1270807), and shared-shift C1/K8 has 27 failed tasks in n=200
(1269529). Their retained logs do not establish the error cause; scheduler
exit codes alone are not enough to attribute these to the reuse guard.

## Interpretation of the three proposed options

The current-DGP `min`/`1se` target-anchor comparison is also complete: 40/40
seeds per configuration. All 120 independently recomputed `min` estimates,
biases and SEs agree with the RCAL experiment's original-anchor rows within
1e-12. Both nuisance models use the same data and CV seed under the two rules.

| Config | min bias | 1se bias | min RMSE | 1se RMSE |
|---|---:|---:|---:|---:|
| C1 | -0.002302 | -0.010196 | 0.033326 | 0.033161 |
| C2 | -0.005572 | -0.012997 | 0.035688 | 0.035309 |
| C3 | -0.007560 | -0.020040 | 0.028765 | 0.033427 |

For C3 the paired estimate shift is -0.01248 (MCSE 0.00113), and the paired MSE
increase is 0.000290 (MCSE 0.000105). The 40-run coverage rises from 0.950 to
0.975 despite the worse bias and RMSE, illustrating why coverage alone is not
an acceptance criterion. These findings concern the target-anchor `1se` change;
they do not evaluate every possible source or aggregation tuning strategy.

Removing all squared terms makes the currently correct C3 outcome working model
misspecified: nonzero quadratic coefficients occur in the true mechanism. It is
not a neutral dimensionality reduction. RCAL is a well-motivated controlled
experiment, but should include arm-specific calibrated PS and matched weighted
OR; its effect on full TATE must be measured after rebuilding aggregation.
Changing the lambda rule is a sensitivity analysis, not a way to select a method
by observed coverage. Earlier one-SE results used an older DGP and cannot be
silently carried over to the current X-dagger design.

For a logistic outcome, the AIPW derivative with respect to its coefficients is
`(1 - I[A=a]/p_a) * m_a * (1-m_a) * phi(X)`. Ordinary calibrated propensity
estimation balances `phi(X)`, which does not by itself balance this derivative's
extra `m_a*(1-m_a)` factor when the propensity model is misspecified. This explains
why standard RCAL is not an automatic C3 inference repair. Tan's model-assisted
interval guarantee allows outcome misspecification when PS is correct; the
additional doubly robust interval statement is for a linear outcome model.
A future target-calibration candidate should explicitly match the relevant
score derivatives and be checked under the actual binary-outcome design.

Primary references: [RCAL manual](https://cran.r-universe.dev/RCAL/doc/manual.html)
and [Tan, model-assisted inference](https://arxiv.org/abs/1801.09817).

## Repair and validation priorities

1. Restore the claimed inner-fold isolation for tuning and propagate the cache
   setting through saved fits and rho reuse. Keep the uncached computation as
   the reference; cache only identical training problems. The observed effect
   on production bias/coverage still needs a paired experiment after that repair.
2. Route final outcome fitting through the shared solver dispatcher and record
   the actual solver. Independently derive baseline inference for the nuisance
   fits actually used; the MLE sign correction alone does not validate lasso SEs.
3. Investigate residual borrowing in shared shift using prespecified alternatives
   such as the existing hard-threshold sensitivity or arm-aware protection of
   common weights. Include the corresponding covariance and weight-learning
   terms. The armwise diagnostic is evidence for this investigation, not an
   automatically approved replacement primary estimator.
4. If the n=1000 nuisance bias remains, test a calibration or regularization-bias
   remedy matched to the binary outcome score, including source nuisances.
   Preserve raw p=100 and the true quadratic signal terms; do not change the DGP
   merely to improve coverage. Target-only RCAL did not resolve this pilot.
5. Freeze a new candidate and evaluate all scenario cells with complete planned
   replicates, paired MSE uncertainty, coverage uncertainty, SE/SD and failures.
   Use independent confirmation seeds after selecting a candidate. Current v4
   failures establish that the present implementation is not uniformly adequate;
   they do not establish that satisfactory performance at n=1000 is impossible.

Full results and logs are indexed by `run_ledger.json` in the scratch root.
The production estimator has not been changed or republished by this audit.
