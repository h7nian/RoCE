# Merged-weight CV: recovery of an exhausted lambda grid

The first W_overlap RHC-model replicate failed at source s2, treated arm,
training folds 1–8. Its original maximum penalty was 0.1133643. The slope
gradient bounds at the intercept-only optima of the five CV training problems
were 0.1657309, 0.1282918, 0.1331595, 0.1487665 and 0.1420216. PN and CD both
failed, with and without the CV certificate. A larger upper endpoint produced
a converged candidate. The exact inputs and all four failures are retained in
`real_data/rhc/model_validation_v1/diagnostics/initial_cv_failure_v1/`.

The issue is that lambda_max for the complete training sample need not cover
every CV training subset. Sparse indicators can disappear from a subset. In
an exponential-tilting objective this may create an unbounded direction at
insufficient slope penalties; a larger iteration limit does not fix that.

Write the initial smooth loss as

\[
  L(\gamma)=\nu^\top\gamma+
       c\,\overline{\exp(-\gamma_0-X^\top\gamma_{-0})},
\]

where c is the source-arm fraction used by the CV normalization, and nu is
the supplied target moment. The intercept-only optimum is
gamma_0=log(c/nu_0). Its slope gradient in a CV training subset is

\[
  G_j=\nu_j-\nu_0\bar X_{j,\mathrm{train}}.
\]

Every subset mean lies between that feature's observed extrema. Therefore

\[
  B=\max_j\max\{
       |\nu_j-\nu_0\min_i X_{ij}|,
       |\nu_j-\nu_0\max_i X_{ij}|\}
\]

bounds all such slope gradients. For a feasible intercept-only optimum,
lambda>B satisfies its slope KKT inequalities. With clipped tilting losses,
this argument also requires the intercept solution to lie inside the clipping
range. The implementation still requires actual convergence and a finite
validation loss in every fold; this bound is not a substitute for those checks.

The recovery rule applies only after the specific error that the entire
original initial-weight CV grid failed. It appends geometrically increasing
penalties up to B times (1+1e-6), retains every original candidate, restores the
pre-CV RNG state, and runs the same CV rule again. It makes only one retry.
Other errors and a failed retry propagate normally. The final fit must pass
the existing convergence checks. The original and extended grids are saved in
the coefficient attribute `cv_grid_retry`.

Successful original CV paths and fixed-lambda fits do not change. No tolerance,
compiled solver, loss, aggregation cutoff, fold boundary or communication round
changes. In this fixture the conservative extrema bound is much larger than
the actual fold-specific bound; several large penalties yield identical
intercept-only fits. The inherited minimum-CV-loss tie rule can select the
largest of these equivalent penalties. This is separate from the aggregation
penalty and from the source truncation radius.

Validation in `implementation/r12/initial_dr_cv_retry_v1/` includes:

- A rare-indicator fixture where the original grid fails, exact preservation
  of the RNG/folds, retention of the old grid, and PN/CD comparisons.
- Successful-path and fixed-penalty checks, plus the existing CV certificate,
  KKT, source-rule and RHC-interface tests: 18 tests and 292 assertions pass.
- Recovery of the exact failing RHC-model subset, with strict final convergence.
- Complete reruns of the real Private/min RHC analysis at source radii 3 and 5:
  estimates, standard errors, nuisance models and aggregation weights reproduce
  their saved references within 1e-10.

The earlier failed study stays frozen. The repaired study uses the same
population template, replication seeds and radii. Its inference performance
must be assessed separately from this numerical recovery check.

The expanded batch subsequently exposed the same defect in calibrated weights
(W_overlap, repeat12, source s3, treated arm, outer fold4). The original upper
endpoint is0.03905091, but CV fold3 has a constant-feature support penalty floor
of0.04120191. Thus every original penalty is below a necessary coercivity bound
in that CV fold. This is stronger evidence than a slow or unconverged optimizer.
The original failure, exact fit arguments, fold bounds and successful expanded
CV result are retained in `real_data/rhc/model_validation_v2/diagnostics/`.

For this derivative-weighted loss, let h_i>=0 denote the initial outcome
derivatives. The intercept-only slope gradient replaces the ordinary mean by
sum(h_i X_i)/sum(h_i). This weighted mean has the same extrema bound B, provided
the derivative sum is positive and the intercept optimum is feasible. The
shared `.density_ratio_cv_with_retry` helper therefore handles initial and
calibrated/refined weights without maintaining duplicate recovery logic.
The new implementation still requires separate regression, exact-fixture and
successful-path parity checks before replacing any saved numerical result.
