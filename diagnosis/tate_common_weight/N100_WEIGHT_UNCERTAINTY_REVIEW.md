# C1/K2 weight uncertainty review

## Latest checkpoint

The final n=100 prefix is complete and independently hash/statistically reviewed
(job 18384884). All 7,200 numerical checks pass; no recorded hard nuisance or
weight-bootstrap failures occurred. Original coverage is 82% at rho=1 and 83%
at rho=1.5; the relearned diagnostic covers 88% at rho=1, with exact pointwise
95% MC upper bound 0.936431. Mean analytic SE/SD is 0.768/0.754 at those rhos.
At rho=1, relearned SE/SD is 0.968 but mean error is +0.010575 (MCSE 0.003047);
at rho=1.5 mean error is -0.000423 (MCSE 0.003332). Target-only itself has mean
error -0.011766 (MCSE 0.003288), coverage 92%, and RMSE 0.034767. All six RoCE
RMSEs remain numerically lower. These are descriptive paired/pointwise findings,
not adjusted sequential tests or validated inference. Saved-data n=100 tail
review 18385881 now passes exact reconstruction of all 600 fits and 1,800 site
cells; high-leverage observations remain present and retained. The final
inference argument remains open; formal production is not
submitted. See the leading table in `WEIGHT_BOOTSTRAP_SIMULATION.md`.
The preceding-stage history below is retained as history, not live-job status.

Saved-data mechanism job 18410137 subsequently verifies the exact bias
decomposition: at rho=1, target -0.011766 plus s1 +0.022767 and s2 -0.000427
gives RoCE +0.010575. The mean s1 fold weight is still 0.14934. At rho=1.5,
subtracting the observed MC mean bias in a truth-assisted thought experiment
leaves analytic coverage at 83%; using the observed MC SD instead gives 95%.
These are descriptive diagnostics, not deployable intervals or new validation.
See `N100_COVERAGE_MECHANISM_DIAGNOSTIC.md` for all 600 reconstruction checks,
6,000 fold/source contributions, paired variance decompositions and limitations.

The later stepwise audit identified actual duplicate-origin overlap in the
row-resampled nuisance CV. A grouped-CV implementation now passes the complete
v19 software gates. The corrected runner now passes both production-dimensional
identity gates and its first actual rho=0 draw's numerical, origin-partition,
and optimizer audits; the first actual rho=1 draw also passes those audits.
Both predeclared 20-draw batches are now complete and audited, while independent-
sampling statistical validation remains pending. The next fresh-seed gate is
predeclared in `INDEPENDENT_INFERENCE_PILOT.md`. Its runner passed 122 targeted
assertions and stubbed launch checks. Task 1 (fresh seed 10001, `18279345_1`)
completed and passed the real six-rho execution/numerical gate and independent
review. Tasks 2--5 also completed and passed the exact n=5 integrity and
statistical recalculation review. The subsequent n=10 checkpoint is also
complete and independently reviewed. The n=25 and now n=50 checkpoints also
pass independent integrity/statistical/tail reviews. At n=50 all 3,600 numerical
checks pass, and RoCE RMSE is numerically below target-only at every rho.
However, analytic coverage is 43/50 at rho=0.5 and 42/50 at rho=1; their exact
pointwise MC upper limits are 0.941808 and 0.928299. Relearned SE/empirical-SD
ratios improve to 0.988--1.036 without validating coverage, which remains
90--94%. Target-only itself has mean error -0.014942 (MCSE 0.004293).
These are not validated inference or prespecified sequential/multiplicity-
adjusted tests. All adverse seeds remain included. The final declared pilot
prefix, tasks 51--100, is submitted as `18353074`; audit `18353177` and summary
`18353219` form its success-dependent chain. Concurrency remains 16; 50 CPUs
per job permit five positive-rho workers, with unchanged five CV threads,
source workers, methods and intervals. Formal production is not submitted. See
`WEIGHT_BOOTSTRAP_SIMULATION.md`,
`GROUPED_CV_IMPLEMENTATION.md` and
`STEPWISE_FUNCTION_REVIEW.md`; the older calibration results below are retained
as versioned diagnostic evidence, not relabeled as validated inference.

All ten predeclared v18 calibration seeds have completed and passed bundle
audits. The canonical ten-seed summary is
`results/direct_tate_mc500_b5000/weight_bootstrap_calibration_v18_summary_n010_audit_v1/`.
This is an implementation/descriptive checkpoint, not an accepted primary
variance policy or a passed coverage gate. These later grouped-CV numerical
gates do not change that interpretation. Full-goal requirements are
tracked in `COMPLETION_CHECKLIST.md`.

## Wald penalty scale and numerical-safeguard regime

The frozen implementation minimizes

`N_all * Var_hat(eta) + sum_j (lambda * t_hat_j - 1)_+ * abs(eta_j)`,

where `t_hat_j = abs(source_hat_j - target_hat) /
sqrt(max(var_discrepancy_j, 1e-6))`. The discrepancy variance includes both
target/source components and their target-side covariance. See
`R/model_fitting.R`'s `optimize_weights` contract,
`R/cross_fitting_aggregation.R`'s Phase-2 Wald calculation and
`src/weight_optimization.cpp`'s penalty construction. The C++ optimizer uses
the unscaled variance objective and divides the main-scale L1 penalty by
`N_all`; this is algebraically equivalent, not a missing sample-size factor.

Two regimes must be distinguished. In the intended statistical regime where
the variance safeguards are inactive, assume compatible-source discrepancy
zero, a joint root-n CLT, consistent discrepancy SEs, fixed site proportions
and a fixed number of folds/sites. Then the standardized null discrepancy
can have a nondegenerate `abs(Z)` limit. A fixed-cutoff penalty
`(lambda * abs(Z) - 1)_+` is consequently random, not necessarily vanishing.
The limiting variance quadratic can be deterministic while its penalized
minimizer remains random. An exception requires additional structure making
the minimizer insensitive to the random penalty; it is not implied by
convexity or numerical convergence alone. A formal asymptotic version of this
regime must specify safeguards tending to zero fast enough (for example a
variance floor `epsilon_n = o(1/n)`), rather than silently retaining a fixed
machine safeguard while assuming an always-consistent SE.

In contrast, the literal code holds `VARIANCE_MIN = VAR_MIN = 1e-6` fixed
(`R/constants.R` and `src/numerical_constants.h`). At mathematical n tending
to infinity, an O(1/n) null discrepancy variance eventually falls below that
floor. Its denominator then stays 0.001 and its null Wald statistic tends to
zero, not `abs(Z)`. Fixed curvature floors and stopping tolerances can also
alter the optimizer's extreme-n limit. That different numerical-algorithm
limit is not a justification for the intended finite-sample Wald procedure.
The initial review's unqualified nondegenerate limit for the literal
fixed-floor algorithm was therefore withdrawn and replaced by this distinction.

A read-only n=25 saved-intermediate audit examined 1,500 Phase-1b discrepancy
variances (25 seeds, six rhos, five outer folds, two sources). All were finite;
the minimum was 0.001881532, over 1,881 times the 1e-6 floor. None triggered
that floor. Their natural diagonal curvatures, twice those variances, likewise
exceed the curvature floor. All 3,750 V_ot/V_t/V_s components exceed 1e-6 and
all 750 fold PSD-ridge values are zero. Thus the current pilot is in the
safeguard-inactive working regime; its observed variability is not explained
by these numerical floors. This check is not another coverage summary or a
new simulation checkpoint.

For the exact Phase-3 estimator, write the source's fold contribution as
`L_kj = p_target,k * (M_kj - T_k) + p_source,jk * D_kj`, using the actual
target/source fold fractions. Relative to fixed reference weights eta0, the
weight-replacement remainder is exactly
`sum_kj (eta_hat_kj - eta0_kj) * L_kj`. Under compatible-source root-n
conditions, `L_kj = O_p(n^-1/2)`. If weight deviations remain O_p(1), this
remainder need not be negligible at the root-n scale. Do not replace it with
mean reported weights times pooled source estimates.

The verified site-centered, realized-weight pseudo-value variance identity
does not by itself establish how that remainder and cross-fold dependence
behave under repeated sampling. Each evaluation fold is excluded from its own
weight fit, but observations contribute to other folds' overlapping training
sets; conditioning on all weights does not automatically preserve independent
evaluation noise. A suitable weight-stability result or random-weight
studentization/CLT argument is still required. Absolute values, positive parts
and L1 kinks also make an ordinary smooth weight influence function, or an
independent additive variance correction, unjustified without extra conditions.
These observations do not prove every current CI invalid, nor validate the
weight-relearning bootstrap. No floor, cutoff, estimator or interval policy
was changed in response to this review.

## Deterministic one-source overlap counterexample

`check_wald_weight_limit.R` checks a specific Wald-soft example against the
installed v19 optimizer, without nuisance fitting or Monte Carlo draws. Set
V_ot=0.5, V_t=V_s=C_ot=0.25 and n_t=n_s=1000 in one optimizer call. Then
`N_all * Var_hat(eta) = Q(eta) = 1 - eta + eta^2`, and at lambda=1 its penalized
solution is exactly

`g(u) = max(0, 0.5 - 0.5 * (abs(u) - 1)_+)`.

The script matches this formula to 25 deterministic optimizer inputs within
1.12e-16. This is a feasible one-coordinate variance quadratic, not a parameter
fit to the n=50 results. The numerical variance floor is inactive in the check.

For a deliberately idealized Gaussian fold model, let independent Z_k and E_k
be standard normal, A_k=-0.5 Z_k+sqrt(0.75) E_k, and
`U_-k = sum_{l!=k} Z_l / sqrt(m-1)`. Define

`L_m = sum_k {A_k + g(U_-k) Z_k} / sqrt(m)`.

The realized-weight proxy is the fold average of `1-g(U_-k)+g(U_-k)^2`.
Its expectation is `1-E[g]+E[g^2] = 0.776294825698`. For distinct folds,
conditioning on their shared complement noise
`R ~ N(0,m-2)` gives cross covariance `c_m=E[q_m(R)^2]`, where, for
`a=sqrt(m-1)`,

`q_m(r) = -{Phi(2a-r)-Phi(a-r)-Phi(-a-r)+Phi(-2a-r)} / (2a)`.

Direct piecewise integration independently checks this integration-by-parts
formula within 3.01e-17. Consequently
`Var(L_m) = E[V_fix] + (m-1)*c_m`:

| Fold count m | E[V_fix] | Pair cross covariance c_m | Var(L_m) |
| ---: | ---: | ---: | ---: |
| 2 | 0.776294825698 | 0 | 0.776294825698 |
| 3 | 0.776294825698 | 0.009135695103 | 0.794566215904 |
| 5 | 0.776294825698 | 0.007419893561 | 0.805974399942 |

The m=2 row is a mathematical comparison, not a supported two-level RoCE
fitting configuration. With m>=3, the shared complement noise prevents the
special even-symmetry cancellation of the two-fold example. The extra variance
is about 2.35%/3.82% of the proxy in these chosen examples. Those numbers are
not estimates of the real pilot's missing variance, an SE inflation rule, or
a proof of the complete RoCE limiting distribution. Matching an average
variance also does not establish normality after random studentization.

## Small-information bootstrap: analysis only

One possible way to address the random empirical drift can be examined in
that abstract moment model. If the original root-n error has representation
`F(G_n)+o_p(1)`, an ordinary nonregular plug-in bootstrap may instead produce
the conditional law of `F(G_n+G_n*)-F(G_n)`, which depends on the realized
training noise. For information size b with r=b/n tending to zero and b
tending to infinity, the corresponding candidate expansion is

`sqrt(b)*(theta_b* - theta_hat_n)
 = F(sqrt(r)*G_n + G_b*) - sqrt(r)*F(G_n) + o_P*(1)`.

Under a suitable conditional joint CLT, tightness, almost-sure continuity of
F and uniform control of the approximation remainder, this tends to `F(G*)`.
F need not be homogeneous: the second term comes from rescaling the original
estimator. Variance convergence additionally needs conditional uniform
integrability of squared errors. The n-scale variance would be
`r * Var*(theta_b* - theta_hat_n)`, and an interval would need justified
conditional quantiles rather than automatically multiplying an SE by 1.96.

Independent positive Gamma(r, scale=1/r) observation weights have mean 1 and
variance 1/r. For centered score vectors h_i, the conditional covariance of
`sqrt(b)/n * sum_i (W_i-1)*h_i` is exactly the empirical covariance
`sum_i h_i h_i^T/n`. Sharing each observation's weight across its fold/arm
appearances retains overlap terms. This is a triangular-array argument as
r tends to zero; the squared standardized individual multipliers are not
uniformly integrable, so a fixed multiplier-distribution theorem cannot simply
be quoted. Lindeberg/no-dominant-observation conditions, consistent empirical
covariances, stable nondegenerate weight objectives and effective sample sizes
for each actual site/training subset would need to be checked. Perturbing only
sample-size denominators without score-mean perturbations supplies no sampling
law. Covariance blocks may be consistently reweighted or consistently held
fixed in an explicitly justified limiting approximation; reweighting them is
not, by itself, a proof of the required CLT.

Actual RoCE also treatment-stratifies its outer folds. An extension must
justify that joint covariance/conditioning structure, not silently identify
it with independent Gaussian folds. The high-dimensional nuisance remainder,
fixed numerical safeguards and finite-sample leverage remain additional
requirements. No b sequence, Gamma option, changed CI or new simulation
protocol has been implemented or selected from these calculations.

### Treatment-stratified folds: what remains to justify

`partition_into_folds` in `R/cross_fitting_algorithms.R` assigns treated and
control observations separately. Fold treatment counts therefore share the
random site treatment total. A simple statistic such as each fold's mean A
has a common treatment-composition fluctuation under unconditional sampling;
independent multipliers on separately normalized, disjoint fold means do not
automatically reproduce that joint distribution. This is a caution about a
candidate resampling argument, not an established error in RoCE's variance.

In particular, the actual `.weight_bootstrap_outer_estimate` does **not** sum
separately normalized fold means. It assembles fold-weighted pseudo-values by
original observation ID and takes one normalized weighted mean per entire
site (lines 533--583 of `R/aggregation_weight_bootstrap.R`). The shared site
denominator matters. Do not apply a zero-crosscovariance argument for disjoint
fold means directly to these globally normalized contributions. Even in a
simple deterministic-coefficient example, fold-marginal and cross-fold errors
can cancel for equal coefficients; random adaptive coefficients need their
own joint-law argument.

Marginal centering also does not imply that score means conditional on A
are zero. But a nonzero empirical arm mean from one saved seed does not prove
a nonzero population conditional mean or a first-order bootstrap error.
There are important possible cancellations: if target and source-assisted
regressions share the same correct target response function, their target-side
discrepancy removes that regression term, leaving residual terms with zero
conditional mean given (X,A). The common target regression term in Phase 3
also has coefficient `(1-sum eta)+sum eta=1`. These are sufficient-structure
examples, not verified assumptions for all robustness regimes in this project.
Outcome misspecification, calibration corrections and site shifts must be
handled explicitly before deciding whether a composition term survives.

A conditional-on-A analysis would need its own influence-function and nuisance
normalization argument, including the distributions of X within treatment
groups. Simply subtracting empirical arm means or normalizing Gamma weights
within arms is not a validated fix. No fold rule, multiplier, point estimator
or CI has been changed based on this concern.

## Local manuscript provenance and rate conditions

Both `docs/main.tex` and the nested repository's `overleaf/main.tex` have
uncommitted modifications. The earlier review's description of the latter as
git-clean/immutable was incorrect and is withdrawn. `JASA/main.tex` is absent
from this checkout. No remote fetch was performed and no manuscript was edited
during this review; local text cannot establish authorship or the latest
remote comments.

The nested Overleaf repository's committed main.tex revision
`f6a9cf2a9682a206877fb897efd772f79c0da706` (2026-07-06) already requires
`lambda_N -> 0` and `lambda_N*sqrt(N_all) -> infinity` for its oracle theorem.
That committed version describes armwise aggregation, a compact weight ball
and a finite cutoff 2; it must not be relabeled as the current common-TATE,
unconstrained, cutoff-1 implementation. The two current working copies state
the same threshold rate and explicitly distinguish their fixed-cutoff
implementation from the rate-indexed theorem. Thus neither version justifies
automatically applying an oracle normal limit to fixed-cutoff pilot intervals.
This is local method-text evidence, not a claim about Hou's original or remote
manuscript instructions.

## Evidence and interpretation

The completed v15 experiment is stored in
`results/direct_tate_mc500_b5000/roce_tate_hard_screen_v2/`.
Its `C1_K2_n100_passed.txt` gate establishes implementation QC for six rho
settings and 100 replications each. It does not establish satisfactory
coverage or authorize the new bootstrap implementation.

| rho | Soft RMSE | Target RMSE | Soft coverage | Soft mean SE / empirical SD | Bias-centered soft coverage |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | 0.022393 | 0.033869 | 0.96 | 1.117 | 0.96 |
| 0.5 | 0.026334 | 0.033869 | 0.91 | 1.031 | 0.94 |
| 1 | 0.033401 | 0.033869 | 0.85 | 0.800 | 0.88 |
| 1.5 | 0.032800 | 0.033869 | 0.88 | 0.786 | 0.86 |
| 2 | 0.030871 | 0.033869 | 0.88 | 0.854 | 0.90 |
| 2.5 | 0.029658 | 0.033869 | 0.89 | 0.907 | 0.91 |

Bias-centered coverage subtracts the mean Monte Carlo error before checking
the reported interval width. It is a diagnostic using simulation truth, not
an available bias correction for real data. Soft RMSE is numerically lower
than target-only in these 100 replications; this is not a guarantee of
population RMSE dominance. At rho 1, the paired MSE difference is only
-0.00003149 with MCSE 0.00019327.

The fixed-weight treated-minus-control pseudo-value construction retains
cross-arm covariance and uses within-site centering. For random learned
weights, the derivative of the aggregation functional with respect to a
source weight is the source-target discrepancy. This creates a possible
additional weight-learning contribution when that discrepancy is nonzero.
Outer sample splitting alone does not account for every observation's
influence through other folds' learned weights. The empirical results are
consistent with this mechanism, but do not uniquely prove it or exclude
nuisance estimation and regularization bias.

## Bootstrap validation status

`estimate_tate_weight_bootstrap()` is an opt-in diagnostic with fixed
nuisance fits and shared observation multipliers. The original analytic SE
remains the reported SE. Hard screening is nonregular near its threshold;
relearning the screen in bootstrap draws is not itself a validity proof.

An end-to-end small-data check at sim_id 901, n_total 600, K 2, p 4,
three folds, RoCE DGP, and three nuisance lambdas produced analytic SE
0.06324911. With 1,000 fixed-weight draws and the simulation's bootstrap
seed for sim_id 901, bootstrap SE was 0.06448059 (ratio 1.01947).
This repeats the same dataset as the earlier 10-draw check; the separate
sim_id 902 check cannot establish Monte Carlo error for sim_id 901.

The v15 n=100 output directory contains CSV summaries, not the fitted raw
inner/outer component RDS objects required by the new bootstrap. Those
summaries cannot support an exact offline weight-relearning bootstrap.
Production-scale calibration must refit the required datasets and retain
the fitted components in a new versioned output directory.

## Follow-up implementation review

The initial bootstrap tests did not cover several important cases. The
follow-up review reproduced and repaired these issues:

- Inner training IDs could incorrectly refer to the outer evaluation fold.
  Validation now checks the exact observation partition and ordered inner
  fold correspondence before generating multipliers.
- New bootstrap failure statuses did not trigger Slurm implementation
  failure. Both malformed diagnostics and failed draws now do so.
- Bootstrap pooling used raw values across inner folds, adding between-fold
  mean variation absent from the production weight objective. It now pools
  the reweighted components centered within inner folds. On the real sim_id
  902 fit, all-one multipliers reproduce variance moments to 2.3e-16 and
  weights to 1.2e-16.
- Incomplete schemas and nonlogical flags could escape diagnostic validation;
  these are rejected. Two-round comparison rows correctly retain missing
  bootstrap diagnostics when one-round bootstrap is enabled.
- Inner-CV selection is explicitly unsupported by this diagnostic until its
  pooling also propagates multipliers. Fixed numeric lambda is supported.
- Generated Rd files now include the new simulation arguments.

After these repairs, the core bootstrap tests have 101 passing assertions;
simulation diagnostics and Slurm policy tests pass, and Rd parsing reports
no warnings. The queued v16 package checks were cancelled before they ran;
v17 subsequently passed both installed-package and R CMD check gates:

- Installed tests: job `18137228`, completed successfully.
- R CMD check: job `18137248`, completed with exact `Status: OK`.
- Source fingerprint:
  `49a1e8e49577328a0e9e2273db504bb69a01e40623e2c725a612c671035f5c78`.
- Installed-package fingerprint:
  `34124883e12ff59296797aa811e653b3df02fc6ef0bc6985c93eaca3ba80dc58`.

Both fingerprints were independently recomputed and matched their gates.
The seed-1 calibration pilot `18138523` subsequently failed after 01:04:32
with `calibration did not retain complete artifacts for all six rhos`.
No result directory was published and the task lock was removed. Its logs
remain under `results/direct_tate_mc500_b5000/logs/18138523_*`.

The failure was traced to group artifact keys: vector `format()` allocated
`0.0`, `1.0`, and `2.0`, while scalar writes used `0`, `1`, and `2`.
Preallocated NULL slots survived alongside appended fitted objects. The
group helper and calibration validator now use scalar-formatted keys
consistently. A regression test covers all six rhos, an out-of-order subset,
and no requested artifacts. This is an artifact handoff correction; it
does not change the TATE point estimator or bootstrap calculation. A new
source version must pass software gates and rerun the pilot before any
production-dimension calibration result can be reported.

## Successful v18 seed-1 calibration

Jobs `18149939` (installed tests) and `18149998` (R CMD check, `Status: OK`)
passed for source fingerprint
`2f6fcf63e0b2143bd8429763d1f1afc4c1b012aadd57c4a88eae7496900c8375`
and installed-package fingerprint
`27ebfb159e6e7f6db955be810e3ef89724741a6f02f68e83730a8f808e774b45`.
Pilot `18150037` completed with exit code 0 in 01:03:51. Its output is
`results/direct_tate_mc500_b5000/weight_bootstrap_calibration_v18/seed_000001`.

The audit verified all four payload hashes, 84 rows (14 methods at six rho
values), zero implementation/bootstrap failures, and six non-null fitted
objects with valid observation-ID partitions. All 84 rows agreed with the
v15 seed-1 outputs across 114 statistical fields to tolerance 1e-12. The
saved fitted estimates and analytic SEs also matched the CSV rows.

| rho | Analytic SE | Fixed-weight bootstrap SE | Relearned-weight SE | Relearned / fixed |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 0.021691 | 0.022352 | 0.023099 | 1.033 |
| 0.5 | 0.022678 | 0.022721 | 0.026601 | 1.171 |
| 1 | 0.025399 | 0.025014 | 0.032540 | 1.301 |
| 1.5 | 0.026439 | 0.026096 | 0.029231 | 1.120 |
| 2 | 0.026439 | 0.026096 | 0.026997 | 1.035 |
| 2.5 | 0.026439 | 0.026096 | 0.026466 | 1.014 |

Fixed-weight bootstrap/analytic SE ratios range from 0.985 to 1.030 at
B=500. The larger relearned SE near rho=1 is consistent with additional
weight-learning variability, but one seed cannot establish sampling coverage
or validity under nuisance misspecification. The primary analytic intervals
remain unchanged. Calibration seeds 2 and 3 were submitted as jobs `18160222`
and `18160232`, respectively, using the same frozen package and configuration.

Replaying the saved rho=1 fit with its recorded seed and B=500 reproduced
both stored bootstrap SEs to tolerance 1e-12. Writing F for the fixed-weight
bootstrap estimate and D for the relearned-minus-fixed estimate gives:

| Bootstrap variance component | Value |
| --- | ---: |
| Var(F) | 0.000625715542 |
| Var(D) | 0.000182096974 |
| 2 Cov(F,D) | 0.000251034750 |
| Var(F+D) | 0.001058847266 |

The variance decomposition identity held to machine precision. The positive
covariance contribution shows why adding a separate weight variance while
ignoring covariance would not reproduce this bootstrap distribution. This
decomposition is conditional on the saved nuisance fits and is not a proof
of unconditional frequentist coverage.

## Seeds 2 and 3

Jobs `18160222` and `18160232` completed successfully in 01:06:27 and
01:07:11. Both output bundles passed all four payload hash checks and had
84 rows, six complete saved fits, zero implementation failures, and zero
bootstrap failures. All 114 checked statistical fields matched their v15
counterparts to tolerance 1e-12. Saved fit estimates and analytic SEs matched
their respective CSV rows.

| rho | Seed 2 relearned / fixed SE | Seed 3 relearned / fixed SE |
| ---: | ---: | ---: |
| 0 | 0.997 | 1.018 |
| 0.5 | 1.089 | 1.260 |
| 1 | 1.341 | 1.301 |
| 1.5 | 1.594 | 1.182 |
| 2 | 1.737 | 1.073 |
| 2.5 | 1.901 | 1.014 |

Seed 2's large high-rho increase needs further statistical review; successful
execution does not validate its inferential interpretation. The three
replications are insufficient for a coverage conclusion.

A same-seed B=500 replay of seed 2 at rho=2.5 reproduced the recorded fixed
SE 0.033391749959, relearned SE 0.063475864922, and ratio 1.900944544666 to
within 3e-15. Its fitted s1 fold weights were approximately
`0, 0, 0.311, 0.403, 0.249`, with Wald statistics
`6.40, 6.39, 3.81, 3.28, 4.20`. Seeds 1 and 3 gave s1 zero weight in all
five folds at this rho. In the bootstrap, seed 2 folds 3--5 sometimes gave
s1 zero weight and sometimes weights above 0.5; all weights remained finite,
nonnegative, and no larger than 0.811. This identifies soft-weight switching
and magnitude variation, without evidence of numerical divergence.

The replay's variance decomposition was approximately 0.0011150 fixed,
0.0013355 relearned-minus-fixed, and 0.0015787 twice their covariance,
summing to 0.0040292. The positive covariance accounts for about 39% of
that bootstrap variance. This is evidence about the fitted bootstrap
distribution, not a validation of its unconditional confidence intervals.

Seed 3's fixed-bootstrap/analytic SE ratio at rho=2.5 was 1.06888 with B=500.
Repeating the fixed-weight bootstrap on the same fit and recorded seed with
B=5,000 gave SE 0.02576500 versus analytic SE 0.02566709 (ratio 1.003814).
The increase in draw count supports Monte Carlo fluctuation as the explanation
for this fixed-weight discrepancy. It does not validate the relearned-weight
distribution or substitute for a full-refit comparison.

Calibration seeds 4 and 5 were submitted as jobs `18170182` and `18170278`
using the same frozen v18 package, six rho settings, and B=500.

Seed 5 completed successfully in 00:52:26. Its four payload hashes passed;
84 rows and six saved fitted objects passed identity and observation-ID
checks, with zero implementation/bootstrap failures. All 114 checked
statistical fields matched v15 to tolerance 1e-12. Relearned/fixed SE ratios
at rho 0, 0.5, 1, 1.5, 2, and 2.5 were respectively
1.035, 1.053, 1.280, 1.292, 1.234, and 1.190. Seed 4 was still running when
this audit was recorded, so no complete five-seed summary was calculated.

Seed 4 subsequently completed successfully in 00:58:02. Its four hashes,
84 rows, six fitted objects and all 114 v15 comparison fields passed the
same checks, with zero implementation/bootstrap failures. The complete
five-seed execution checkpoint therefore passed for all 30 seed-rho
combinations. This is an execution checkpoint, not a coverage gate.

| rho | Median relearned / fixed SE | Minimum | Maximum |
| ---: | ---: | ---: | ---: |
| 0 | 1.018 | 0.997 | 1.035 |
| 0.5 | 1.089 | 1.048 | 1.260 |
| 1 | 1.301 | 1.234 | 1.341 |
| 1.5 | 1.292 | 1.120 | 1.594 |
| 2 | 1.234 | 1.035 | 1.737 |
| 2.5 | 1.187 | 1.014 | 1.901 |

The remaining predeclared calibration seeds 6--10 were queued in two
success-dependent chains, limiting this batch to two concurrent jobs:

- Seed 6: `18177644`; seed 8: `18177646` after seed 6; seed 10: `18177761`
  after seed 8.
- Seed 7: `18177645`; seed 9: `18177647` after seed 7.

All use the same v18 installation, C1/K2, p=100, six rhos and B=500. A failed
parent prevents its successor from starting. Completed bundles still require
hash, provenance and statistical review before any larger experiment.

Seeds 6, 7 and 9 subsequently completed as jobs `18177644`, `18177645` and
`18177647`, in 01:24:01, 01:02:39 and 01:03:21 respectively. Their four payload
hashes, 84-row/six-rho/14-method schemas, six complete artifacts, and fit-to-CSV
estimate/analytic-SE identities passed. Observation-ID validation was also
explicitly rerun for seeds 6 and 9. All three had zero implementation or
bootstrap failures. The independent v15 comparison covered 436 nonempty
invariant common fields at tolerance 1e-12; it did not claim to reuse the
earlier named 114-field audit list. Thus seeds 1--7 and 9 have passed bundle
audits (48 seed-rho combinations); this is not the predeclared ten-seed gate.
At the last authoritative poll, seed 8 was running and seed 10 was waiting
on its success. No partial eight-seed coverage conclusion is drawn.

Seed 8 then completed as job `18177646` in 01:05:13 and passed the same bundle,
observation-ID, fit-CSV, bootstrap/QC and 436-invariant-field v15 checks, with
zero failures. Seeds 1--9 are now audited; seed 10 (`18177761`) started on
the completed dependency and remains the only unfinished calibration seed.

Seed 10 subsequently completed as job `18177761` in 01:00:12. All four payload
hashes, 84-row/six-rho/14-method schema, six complete fitted artifacts,
observation-ID validation, fit-to-CSV identities and zero implementation/
bootstrap-failure checks passed. Its 114 checked v15 statistical invariants
also matched. The summary over the **entire prespecified 1--10 set** passed
its four payload hashes, preserving input and summary-script fingerprints.

| rho | Mean of relearned/fixed SE ratios | Analytic coverage | Relearned diagnostic coverage |
| ---: | ---: | ---: | ---: |
| 0 | 1.046 | 10/10 | 10/10 |
| 0.5 | 1.084 | 9/10 | 9/10 |
| 1 | 1.239 | 8/10 | 9/10 |
| 1.5 | 1.286 | 8/10 | 9/10 |
| 2 | 1.253 | 8/10 | 9/10 |
| 2.5 | 1.219 | 8/10 | 9/10 |

The fixed-weight bootstrap coverage equals analytic coverage at these six
settings. Exact two-sided 95% binomial intervals are approximately [0.6915,1]
for 10/10, [0.5550,0.9975] for 9/10, and [0.4439,0.9748] for 8/10. A plug-in
coverage MCSE of zero at 10/10 must not be interpreted as certainty.

Point estimates are unchanged. In these first ten seeds, soft RMSE is
0.04170 at rho=1 versus target-only 0.03282, and remains above target-only
through rho=2.5. This is not a new point-estimation regression: the same
statistical fields match the archived v15 run. The larger n=100 point summary
above has different, more stable averages. Changing an SE cannot improve
point-estimate RMSE or remove residual borrowing bias. Neither the favorable
SE increase nor this small-set RMSE comparison establishes population dominance.

## Remaining checks before statistical expansion

The read-only `summarize_weight_bootstrap_calibration.R` now produces an
explicit-seed checkpoint with payload checksums, complete schemas, paired
error/coverage summaries, and both input-bundle and summary-code provenance.
The audited five-seed output is
`results/direct_tate_mc500_b5000/weight_bootstrap_calibration_v18_summary_n005_audit_v2/`.
Its four payload hashes passed verification; n=5 coverage/empirical-SD numbers
are descriptive only. Earlier n005 and audit_v1 outputs are retained, not the
current canonical checkpoint.

A separate full nuisance-refit diagnostic has been implemented outside the
frozen package/workflow. It shares site-by-original-fold multinomial counts
between fixed/refitted nuisances and original/relearned TATE weights. Component
and saved-fit checks passed 131 assertions. The first identity refit (seed 1,
rho 0) was submitted as job `18183278`; no full-refit gate has yet passed.
See `FULL_REFIT_REVIEW.md` for the frozen script hashes, execution contract,
and important saved-design, nuisance-CV duplicate-ID, and nonregularity limits.

The exact finite-sample aggregation was separately reconstructed using the
actual target/source fold fractions. All six seed-1 estimates matched to
5.6e-17. The returned `weights` are mean fold weights, not a vector that can
be applied to already-pooled `source_estimates` to reproduce the estimator.
See `AGGREGATION_OUTPUT_CONTRACT.md` for the precise reporting contract,
adaptive-weight decomposition and next-freeze documentation corrections.
Expanded diagnostic/saved-fit tests passed 143 assertions. The package source
hash was rechecked against both existing v18 software gates and remained
`2f6fcf63e0b2143bd8429763d1f1afc4c1b012aadd57c4a88eae7496900c8375`.

The full nuisance-refit identity job `18183278` subsequently completed and
passed independent point/SE/variance/weight/Wald/ID audits exactly (matched
fixed-nuisance weight error at most 5.55e-17). The first actual multinomial
resample and the rho=1 identity are now jobs `18199028` and `18199029`.
Neither is a passed resampling or coverage gate yet. A separate read-only RHC
preflight is recorded in `REAL_DATA_FINAL_PREFLIGHT.md`; it confirms local
data availability but no final variance-frozen RHC result exists.

- Check fixed-weight versus analytic SE and weight-relearning behavior at
  production dimensions, retaining components for reproducible calibration.
- Assess coverage over a predetermined replication set and compare a bounded
  full nuisance-refit resampling reference before adopting a new primary SE.

Full simulation configurations, 500-replication runs, sensitivities, and the
final real-data analysis remain outstanding.
