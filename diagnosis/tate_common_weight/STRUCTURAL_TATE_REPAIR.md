# Structural TATE repair: current priority and acceptance criteria

## User constraint

The user explicitly requires a solution at the estimator/weight-learning/
analytic-variance level. Bootstrap, an SE inflation coefficient, and
truth-assisted recentering are NOT repair candidates. Existing resampling
outputs remain historical diagnostics only. Do not launch new bootstrap
experiments or promote bootstrap CIs as the solution. Operational capture
integration is paused while structural estimation issues take priority.

No root package or manuscript was changed during this structural review.
The frozen v19 installation, all previous data/results, and JASA/Overleaf
content remain untouched. The working Overleaf manuscript is not a pristine
teacher-authored archive and is read only as current local method context.

## 1. Weight stability and residual source bias

The current Phase-2 rule uses `(lambda * Wald - 1)_+` with lambda=1.
The local oracle theorem instead assumes lambda_N -> 0 and
lambda_N * sqrt(N) -> infinity, plus separated incompatible sources and
consistent primitive moments. A fixed cutoff does not satisfy that condition.
The previously checked one-source example shows that null-compatible weights
can remain random on the first-order scale when numerical safeguards are
inactive. Matching the realized-weight variance algebra does not establish
the required asymptotic linear expansion.

The original FACE paper has a rate-indexed adaptive penalty, with a squared-
discrepancy construction and a truncated-Wald alternative
([arXiv v4, Sections 3--4](https://arxiv.org/pdf/2112.09313)). Its conditions
are not a license to replace a rate with a fixed finite-sample cutoff and
retain the same oracle claim. No tuning rate or constant has been selected
from the n=100 coverage results.

A structural replacement must establish that compatible-source weights have
an appropriate stable limit and that incompatible-source contributions are
negligible at the root-N scale under explicitly stated separation conditions.
Local alternatives are a distinct regime; do not promise uniform oracle
selection or automatically valid zero-bias Wald intervals there.

A smooth quadratic bias-weight candidate now has a conditional rate argument
and one fixed exploratory reweighting of all saved C1/K2 n=100 data. It changes
the weight objective, not the SE. Power=3/4 was fixed before evaluation and
has not been tuned. All 600 attempts succeed and point RMSE improves in several
settings, but analytic coverage remains 83--93%; the candidate is not adopted.
See `QUADRATIC_BIAS_WEIGHT_CANDIDATE.md` and
`QUADRATIC_BIAS_WEIGHT_N100_REVIEW.md`. Its full weight-learning derivative and
the existing nuisance/inner-outer alignment gaps remain to be addressed.
The subsequent analytic weight-layer derivative passes actual directional
checks and improves the same candidate's exploratory n=100 coverage from
83% to 93% at rho1 and 86% to 94% at rho1.5, without changing point estimates,
resampling, or applying an SE factor. This calculation still fixes nuisance
fits. A separately derived treatment-balanced design map passes 12 tests and
one saved-seed integration check. See `ANALYTIC_WEIGHT_LAYER_REVIEW.md` for
the derivative, cross-covariance terms, results and remaining inference gaps.

## 2. Calibration equations versus the full nuisance tangent

There is now a concrete, deterministic counterexample to inferring full
outcome-score orthogonality solely from convergence of the implemented
dual-basis calibration equations.

Code references:

- `R/data_generation_face.R:561`: C3 uses linear Z_site but quadratic W_outcome;
  C2 has the opposite basis inclusion pattern.
- `R/model_fitting.R:525`: `.mean_glm_gradient_site_basis` returns
  E_t[g'(W alpha) Z], not E_t[g'(W alpha) W].
- `src/density_ratio.cpp:84`: the native calibrated tilting moment has Z's
  dimension. The use of W to evaluate g' does not add the missing W directions.

For a correct outcome mean m_alpha and source treated mass r, the population
source-assisted mean is

    theta(alpha,gamma) = sum_x q(x) m_alpha(x)
                      + sum_x r(x) exp(-Z(x)'gamma) [m_true(x)-m_alpha(x)].

At alpha_true it equals the target mean for any finite gamma. Nevertheless,
its derivative in the outcome direction is

    sum_x [q(x)-r(x) exp(-Z(x)'gamma)] g'(W(x)'alpha_true) W(x).

The implemented calibration can zero this expression projected onto Z
without zeroing components of W outside span(Z). Pointwise double robustness
at the true outcome function is therefore not enough to discard all
first-order estimation effects from alpha_hat under an incompatible tilting
working model.

`check_calibration_tangent_span.R` verifies this using the actual installed
v19 C++ solver and an independent deterministic optimizer. For three support
points with a correct quadratic logistic outcome and reduced linear tilting:

- implemented calibration moment error: 7.1783e-11;
- native versus independent optimizer parameter difference: 8.3877e-9;
- uncovered quadratic outcome-score derivative: 0.04865675341;
- finite differencing reproduces that derivative to numerical precision.

The example uses no resampling, nuisance CV, random data or SE correction.
All outcome logits are inside M=5, so outcome truncation is inactive.
Script SHA-256: `af80c6b30c62816d3a48a77511b17bc46ae76df297da5084f3ab4ebe78ed95fd`.
This is a C3-style structural counterexample, NOT an explanation of the C1
n=100 shortfall (C1 has matching feature spans). The symmetric C2 direction,
the target anchor, and the actual binary-calibration derivative conditions
must be reviewed separately.

Do not silently replace both nuisance bases by their union and then claim the
original misspecification experiments still test the same robustness regimes.
Any repair must either enforce the needed full tangent equations with a
justified estimator construction or account analytically for the remaining
first-order nuisance effects. Changing working-model definitions is a
method/design change that must be explicit.

## 3. Target anchor and analytic variance

Current target-side evidence is in `TARGET_NUISANCE_SYSTEM_REVIEW.md`: all 75
CV states were retained with zero saved-prediction error, and a real
603-dimensional shared-propensity/two-outcome system passes derivative and
KKT checks. Weak projection regularization still produces a large held-out
shift; no fraction is selected. The full n=100 balanced-design weight-layer
review is also complete, with exploratory coverage 91--94%; see
`BALANCED_WEIGHT_LAYER_N100_REVIEW.md`. Neither result establishes full inference.

A constructive finite-dimensional joint-nuisance score prototype is now
verified against all four installed native fitting steps. It also reproduces
the C2 counterpart and zeroes the complete stacked nuisance derivative in
the three population examples without merging bases. See
`JOINT_NUISANCE_SCORE_PROTOTYPE.md` for the derivation and explicit gaps.
It is not yet a high-dimensional TATE implementation or an aggregation repair.
The subsequent exact joint-observation TATE check also passes both-arm
orthogonality, arm-covariance and unequal-site scaling identities for one
target/source pair. This does not resolve learned multi-source weights or
validate the actual p=100 inference regime.
A regularized blockwise coefficient solver now passes 31 assertions and the
native-kernel joint-observation TATE integration check at zero penalty. The
actual p=100 states are available for a no-refit matrix gate; initial per-k2
states and fold-summed final fits must be distinguished. See
`SPARSE_PROJECTION_IMPLEMENTATION.md`. No production penalty, high-dimensional
inference guarantee or new common-weight rule has been adopted.
The actual p=100 ten-block matrix gate subsequently completed with native
score/KKT and directional-derivative checks; see `SAVED_P100_PROJECTION_GATE.md`.
This is a no-refit implementation gate, not an independent-sampling statistical
gate or evidence that the unresolved C1 coverage problem has been repaired.
The subsequent disjoint-record holdout check shows marked score growth for
the weak-penalty projection despite passing KKT: source SD 0.839 -> 3.022 and
maximum centered source score 4.47 -> 24.89. All probe fractions are retained,
with no holdout-based selection; see `SAVED_PROJECTION_HOLDOUT_REVIEW.md`.
This candidate still needs a principled projection policy and stability
control before any production implementation or statistical expansion.

The cross-fitted target anchor currently fits ordinary penalized propensity
and outcome models. Merely selecting use_rcal=TRUE does not change this path:
`estimate_target_only` documents that it is ignored under use_crossfit=TRUE.
The source calibration theorem cannot simply be assigned to this independent
target estimator. Its nuisance remainder and influence expansion need their
own conditions or an estimator-level repair.

After a coherent estimator/weight construction is established, derive the
TATE influence function with common arm weights, shared-target covariance,
actual site/fold fractions and all nonnegligible nuisance terms. The analytic
SE must follow from that derivation, not from an empirical multiplier chosen
to repair the observed coverage. Then validate the new method with explicit
software and independent statistical gates. The full C1--C3/K2,4,8 primary,
sensitivity and RHC scope in COMPLETION_CHECKLIST.md is unchanged.
