# Baseline Formula-to-Code Mapping (RoCE)

This note documents how each baseline estimator is implemented and where the key formulas appear in code.

## 1) Target-only AIPW

- Estimand: $\mu_t^a = E_t[Y(a)]$.
- Pseudo-outcome:
  - $\phi_i = \hat m_a(X_i) + I(A_i=a)\{Y_i-\hat m_a(X_i)\}/\hat p_a(X_i)$.
- IF variance:
  - $\hat V_{ot} = n^{-1}\sum_i \hat\varphi_i^2$.
  - $\widehat{\mathrm{Var}}(\hat\mu_t^a)=\hat V_{ot}/n$ (with nuisance-adjusted IF in helper).

Code:
- `estimate_target_only`: [R/estimators_target.R](../R/estimators_target.R)
- Cross-fit target-only: [R/estimators_target.R](../R/estimators_target.R)
- AIPW IF helper: [R/estimators_helpers.R](../R/estimators_helpers.R)

## 2) Sample-size weighted (SS)

- Point estimate:
  - $\hat\mu_{SS}=\sum_j (n_j/N)\hat\mu_j$.
- Heterogeneity-adjusted variance (RE):
  - $\widehat{\mathrm{Var}}_{RE}=\sum_j (n_j/N)^2\{\widehat{\mathrm{Var}}_j+\hat\tau^2\}$.

Code:
- Estimator: [R/comparison_methods.R](../R/comparison_methods.R)
- DL heterogeneity (`Q`, $\tau^2$, $I^2$): [R/estimators_helpers.R](../R/estimators_helpers.R)

## 3) Inverse-variance weighted (IVW)

- Point estimate:
  - FE: weights $\propto 1/\hat V_j$.
  - RE: weights $\propto 1/(\hat V_j+\hat\tau^2)$ when $\hat\tau^2>0$.
- Variance:
  - RE variance plus pooled IF safeguard:
  - final variance is `max(var_random_effects, var_pooled)`.

Code:
- Estimator: [R/comparison_methods.R](../R/comparison_methods.R)

## 4) Tilted-AIPW baseline

- Per-site source estimate uses density-ratio weighted AIPW.
- Adds first-order density-ratio score correction (`adj_alpha`) to source IF.
- Final variance decomposition:
  - `target_variance_component + source_variance_component`.

Code:
- Estimator + decomposition: [R/comparison_methods.R](../R/comparison_methods.R)
- Decomposition terms: [R/comparison_methods.R](../R/comparison_methods.R)

## 5) Federated-DR baseline

- Source sites:
  - DR weights via exponential tilting against target moments.
  - Site AIPW IF + first-order DR-weight correction.
  - The density-ratio score/Jacobian uses the full source sample, matching the
    all-observation (`A_dummy = 1`) density-ratio fit; the treatment-arm
    indicator enters through the AIPW pseudo-outcome, not the tilting score.
- Aggregation:
  - IVW point estimate from site-level variances.
  - Variance decomposition consistent with shared-target structure:
    - `var_target_component + var_source_component`.

Code:
- Estimator: [R/comparison_methods.R](../R/comparison_methods.R)
- Shared-target/source decomposition: [R/comparison_methods.R](../R/comparison_methods.R)

## 6) Pooled-DR baseline

- Pools all observations with site-specific DR weights.
- Uses weighted AIPW IF from helper.
- Adds first-order DR-weight correction per source (source + shared-target parts).
- Final variance uses site-stratified within-group centering:
  - $\widehat{\mathrm{Var}} = N^{-2}\sum_g\sum_{i\in g}(\hat\varphi_i-\bar\varphi_g)^2$.

Code:
- Estimator: [R/comparison_methods.R](../R/comparison_methods.R)
- Stacked correction + WSS variance: [R/comparison_methods.R](../R/comparison_methods.R)

## 7) Shared helper pieces

- DR lambda CV for density-ratio model:
  - [R/estimators_helpers.R](../R/estimators_helpers.R)
- DR weights `exp(-eta)` with normalization/clipping:
  - [R/estimators_helpers.R](../R/estimators_helpers.R)
  - Normalized weights are clipped to `[0.1, 10]`; any clipping emits a run-log
    warning so overlap stress is not silent.
- Weighted AIPW now uses nuisance-adjusted IF helper:
  - [R/estimators_helpers.R](../R/estimators_helpers.R)

## 8) Notes on strictness

- Implementations are aligned at first-order IF/asymptotic level.
- Finite-sample exactness is not guaranteed for any semiparametric baseline.
- The code avoids heuristic `dr_inflation` and uses explicit IF-based variance paths.

### Redundant nuisance parameters

The nuisance derivative solver distinguishes numerical redundancy in the active
training design from unstable curvature. It centers and scales the design for
rank checks and transforms sensitivities, data scores and penalty derivatives
in the same coordinates. Full-rank, well-conditioned fits retain the original
direct solve. Redundant fits use the identifiable subspace only when the
sensitivity, score blocks and penalty gradient have negligible components in
the null space. This does not refit the nuisance model or change its selected
lambda, predictors or feature map.

OR checks include the evaluation design, since a relation that holds within one
treatment arm may fail in the evaluation population. Density-ratio checks
include both source and target information. A dependency that fails these
checks, a near-collinear design, or unstable curvature in the retained space
raises an error; no ridge or unverified singular-value truncation is applied.
The numerical rank threshold is `100 * .Machine$double.eps * active_columns`,
the retained relative singular values must be at least `1e-7`, and the null-space
compatibility tolerance is `1e-8`. These are numerical checks, not statistical
tuning parameters or guarantees through active-set transitions.

DR results expose aggregate `curvature_diagnostics` for each arm, site and
nuisance component. The real-data runner writes `curvature_diagnostics.csv` on
success, and `status.txt` and `warnings.txt` even if a method fails. Diagnostics
contain ranks, lambdas and solver checks, not patient-level scores or feature
names. Existing successful outputs are not overwritten: use a new output
directory after reinstalling an updated package. Inference remains conditional
on the selected nuisance lambdas and active sets; Federated-DR also conditions
on its estimated inverse-variance aggregation weights.
Simulation manifests identify this implementation as `conditional_active_set_v2`
so checkpoints using the older derivative solver are not reused silently.
